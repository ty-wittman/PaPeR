#ifndef PAP_NATIVE_LOCUS_EXECUTOR_H
#define PAP_NATIVE_LOCUS_EXECUTOR_H
#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
namespace pap_internal {
struct NativeNoopWorkerGuard {};
struct NativeLocusExecutorJob {
  void* payload = nullptr;
  void (*run)(void*, int) = nullptr;
  const char* label = "locus";
};

template<class WorkerGuard> struct NativeLocusWorkerExecutor {
  const int requested_workers;
  const int effective_workers;
  const int worker_thread_count;

 private:
  std::vector<std::thread> threads;
  std::mutex mutex;
  std::condition_variable job_condition;
  std::condition_variable completion_condition;
  NativeLocusExecutorJob job;
  std::size_t generation = 0;
  int completed = 0;
  bool shutdown = false;
  std::vector<std::string> errors;

  void worker_loop(int worker_id) {
    WorkerGuard worker_guard;
    std::size_t seen_generation = 0;
    while (true) {
      NativeLocusExecutorJob current;
      std::size_t current_generation = 0;
      {
        std::unique_lock<std::mutex> lock(mutex);
        job_condition.wait(lock, [&]() {
          return shutdown || generation != seen_generation;
        });
        if (shutdown) return;
        current = job;
        current_generation = generation;
      }
      std::string error;
      try {
        if (current.run == nullptr)
          throw std::runtime_error("native locus executor received an empty job");
        current.run(current.payload, worker_id);
      } catch (const std::exception& caught) {
        error = caught.what();
      } catch (...) {
        error = "native locus worker raised a non-standard exception";
      }
      {
        std::lock_guard<std::mutex> lock(mutex);
        errors[static_cast<std::size_t>(worker_id)] = error;
        ++completed;
        seen_generation = current_generation;
        if (completed == effective_workers) completion_condition.notify_one();
      }
    }
  }

 public:
  NativeLocusWorkerExecutor(int requested_workers_, int loci)
      : requested_workers(requested_workers_),
        effective_workers(std::min(requested_workers_, std::max(1, loci))),
        worker_thread_count(effective_workers > 1 ? effective_workers : 0),
        errors(static_cast<std::size_t>(effective_workers)) {
    if (requested_workers < 1)
      throw std::invalid_argument("native locus executor requires at least one worker");
    completed = effective_workers;
    if (effective_workers <= 1) return;
    threads.reserve(static_cast<std::size_t>(effective_workers));
    try {
      for (int worker_id = 0; worker_id < effective_workers; ++worker_id)
        threads.emplace_back(
          &NativeLocusWorkerExecutor::worker_loop, this, worker_id);
    } catch (...) {
      stop_threads();
      throw;
    }
  }

  ~NativeLocusWorkerExecutor() { stop_threads(); }

  void stop_threads() noexcept {
    if (threads.empty()) return;
    {
      std::lock_guard<std::mutex> lock(mutex);
      shutdown = true;
    }
    job_condition.notify_all();
    for (std::thread& thread : threads)
      if (thread.joinable()) thread.join();
    threads.clear();
  }

  void dispatch(const NativeLocusExecutorJob& next_job) {
    if (effective_workers <= 1)
      throw std::logic_error("serial locus work must retain its caller-thread path");
    {
      std::lock_guard<std::mutex> lock(mutex);
      if (completed != effective_workers)
        throw std::logic_error("native locus executor dispatched while a job is active");
      job = next_job;
      errors.assign(static_cast<std::size_t>(effective_workers), std::string());
      completed = 0;
      ++generation;
    }
    job_condition.notify_all();
  }

  void wait() {
    if (effective_workers <= 1) return;
    std::unique_lock<std::mutex> lock(mutex);
    completion_condition.wait(lock, [&]() {
      return completed == effective_workers;
    });
    for (int worker_id = 0; worker_id < effective_workers; ++worker_id) {
      const std::string& error = errors[static_cast<std::size_t>(worker_id)];
      if (!error.empty())
        throw std::runtime_error(std::string("native ") + job.label + " worker " +
          std::to_string(worker_id + 1) + " failed: " + error);
    }
  }
};

}
#endif
