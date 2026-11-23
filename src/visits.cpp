#include "visits.h"
#include "decisions.h"
#include "scoring.h"
#include "result_columns.h"
namespace pap_internal::search {

List run_visits(
    SEXP compiled_sexp, List inputs) {
  XPtr<EngineState> pointer(compiled_sexp);
  if (pointer.get() == nullptr) engine_stop("simp compiled snapshot pointer is null");
  EngineState& compiled = *pointer;
  if (!inputs.containsElementNamed("parents"))
    engine_stop("whole-pass input is missing 'parents'");
  const IntegerMatrix parents = inputs["parents"];
  if (parents.nrow() != compiled.n || parents.ncol() != 2)
    engine_stop("whole-pass parent matrix must be n by 2");
  SearchInputs input(inputs, compiled.n);


  input.validate();
  if (input.message_revision != compiled.message_revision)
    engine_stop("whole-pass message revision does not match compiled snapshot");
  SearchState state(compiled, input, parents);
  if (!state.initial_graph_valid())
    engine_stop("whole-pass starting parent matrix violates typed graph constraints");

  // A dense-window benchmark rebuilds the immutable round-start message
  // snapshot, then replays the accepted prefix through this same transactional
  // graph-patch path.  The explicit old state and revisions make accidental
  // replay against a different prefix fail closed.
  for (int row = 0; row < input.replay_moves.nrow(); ++row) {
    const int child = input.replay_moves(row, 0);
    const int old_dam = input.replay_moves(row, 1);
    const int old_sire = input.replay_moves(row, 2);
    const int new_dam = input.replay_moves(row, 3);
    const int new_sire = input.replay_moves(row, 4);
    const int graph_before = input.replay_moves(row, 5);
    const int graph_after = input.replay_moves(row, 6);
    if (child < 1 || child > input.n ||
        state.pass.dam[child] != old_dam ||
        state.pass.sire[child] != old_sire ||
        state.graph_revision != graph_before ||
        graph_after != graph_before + 1)
      engine_stop("whole-pass replay prefix does not match the current graph state");
    if (state.proposal_reason(child, new_dam, new_sire) != 0)
      engine_stop("whole-pass replay prefix contains an infeasible parent patch");
    state.apply_move(child, new_dam, new_sire);
    if (state.graph_revision != graph_after)
      engine_stop("whole-pass replay prefix graph revision did not advance exactly once");
  }

  VisitScorer worker_team(compiled, state, input);
  LocusContext& worker = worker_team.primary();
  std::vector<std::unique_ptr<LocusWorkspace>>& workspaces =
    worker.workspaces;
  EngineCounters& progress_telemetry = worker.pass.telemetry;

  const auto started = std::chrono::steady_clock::now();
  std::vector<AcceptedMove> moves;
  std::vector<VisitSummary> summaries;
  const int focal_end_position = input.focal_end_position == 0 ?
    input.n : input.focal_end_position;
  const std::size_t focal_window_size = input.diagnostic_visit_order.empty() ? static_cast<std::size_t>(
    focal_end_position - input.focal_start_position + 1):input.diagnostic_visit_order.size();
  summaries.reserve(focal_window_size);

  std::vector<ConfigurationDiagnostic> configuration_diagnostics;
  std::vector<ConfigurationLocusDiagnostic>
    configuration_locus_diagnostics;
  std::vector<GraphStateDiagnostic> graph_state_diagnostics;

  std::vector<DiagnosticVisitRow> diagnostic_visit_rows;
  std::vector<DiagnosticPanelRow> diagnostic_panels;
  if(!input.diagnostic_visit_order.empty())diagnostic_visit_rows.reserve(focal_window_size);

  const bool requested_message_model = input.message_scoring;
  for (std::size_t visit_index=0;visit_index<focal_window_size;++visit_index) {
    report_visits(input, visit_index, focal_window_size, moves.size());
    const bool diagnostic_sequence=!input.diagnostic_visit_order.empty();
    const int child=diagnostic_sequence ? input.diagnostic_visit_order[visit_index]:
      input.focal_order[input.focal_start_position-1+visit_index];
    const std::size_t focal_position=diagnostic_sequence ? static_cast<std::size_t>(
      std::find(input.focal_order.begin(),input.focal_order.end(),child)-input.focal_order.begin()):
      static_cast<std::size_t>(input.focal_start_position-1)+visit_index;
    if(diagnostic_sequence && !input.diagnostic_visit_recheck_order.empty()) {
      input.recheck_only=input.diagnostic_visit_recheck_order[visit_index]!=0;
      input.recheck_focals.clear();if(input.recheck_only)input.recheck_focals.push_back(child);
    }
    if (!input.is_recheck_focal(child)) continue;
    state.prepare_visit_ages(child);
    input.birth_lo=state.visit_age_bounds.lo;input.birth_hi=state.visit_age_bounds.hi;
    pap_internal::typed_age::Scope age_scope(&state.visit_age_bounds);



    if(pap_internal::native_demand::live)pap_internal::native_demand::live->begin_visit(child,state.pass.dam,state.pass.sire,state.graph_revision,input.recheck_only,inputs);
    DiagnosticVisitRow* diagnostic_row=nullptr;
    if(diagnostic_sequence) {
      diagnostic_visit_rows.emplace_back();diagnostic_row=&diagnostic_visit_rows.back();
      diagnostic_row->index=visit_index+1;diagnostic_row->focal=child;diagnostic_row->graph_before=state.graph_revision;
    }

    bool message_fallback = false;
    // A gate-window visit that finds nothing assignable while rows waited
    // outside the window is chosen again with every row answered.
    bool window_widened = false;
    std::set<std::pair<int,int>> completed_identity_competitors;
    pap_internal::pair_rescue::CompletionScope completion_scope(completed_identity_competitors);
    int identity_completion_round=0;
    input.gate_window_ln = input.gate_window_base;
    const auto frozen_focal_started = std::chrono::steady_clock::now();
  retry_frozen_focal_model:
    {

    const auto trial_started = std::chrono::steady_clock::now();
    input.message_scoring = requested_message_model && !message_fallback;
    state.pass.message_scoring = input.message_scoring;
    const double message_declines_before =
      progress_telemetry.message_partition_unavailable;
    const bool diagnostic_focal =
      static_cast<int>(focal_position + 1) >= input.diagnostic_start_position &&
      (input.diagnostic_end_position == 0 ||
       static_cast<int>(focal_position + 1) <= input.diagnostic_end_position);
    if (input.configuration_diagnostics && diagnostic_focal &&
        (graph_state_diagnostics.empty() ||
         graph_state_diagnostics.back().graph_revision != state.graph_revision)) {
      GraphStateDiagnostic snapshot;
      snapshot.graph_revision = state.graph_revision;
      snapshot.dam = state.pass.dam;
      snapshot.sire = state.pass.sire;
      graph_state_diagnostics.push_back(std::move(snapshot));
    }
    const auto focal_started = frozen_focal_started;

    state.pass.set_focal(child, state.graph_revision);
    {
      worker.synchronize_focal(state.graph_snapshot(), child);
    }

    const int graph_before = state.graph_revision;
    const int starting_dam = state.pass.dam[child];
    const int starting_sire = state.pass.sire[child];
    const int dam_begin = input.dam_offsets[child - 1];
    const int dam_end = input.dam_offsets[child];
    const int sire_begin = input.sire_offsets[child - 1];
    const int sire_end = input.sire_offsets[child];
    const int dam_count = dam_end - dam_begin;
    const int sire_count = sire_end - sire_begin;
    // Current Sequoia's ordinary Parentage() pass does not reconsider an
    // individual that already has two typed parents.  CheckPedigree() is a
    // separate boundary before the round loop; complete pairs that survive
    // that boundary remain closed during ordinary parentage passes.
    if (starting_dam > 0 && starting_sire > 0 && !input.recheck_only) {

      VisitSummary summary;
      summary.focal = child;
      summary.starting_dam = starting_dam;
      summary.starting_sire = starting_sire;
      summary.selected_dam = starting_dam;
      summary.selected_sire = starting_sire;
      summary.rejection_reason = 5;
      summary.graph_revision_before = graph_before;
      summary.graph_revision_after = graph_before;
      summary.dam_candidate_count = dam_count;
      summary.sire_candidate_count = sire_count;
      summary.complete_pair_skipped = 1;
      summary.message_scoring_model = -1; // No competing likelihood was scored.
      summary.raw_message_approximate = -1;
      summary.focal_wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - focal_started).count();
      summaries.push_back(summary);

      continue;
    }


    const bool incumbent_full_visit = input.incumbent_full_scoring &&
      input.relationship_mode == 1 && (starting_dam > 0 || starting_sire > 0);
    std::vector<ConfigurationScore> configurations;
    const auto active_configuration_started =
      (input.benchmark_opt5_diag_mode || input.bounded_kernel_profile) ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point();
    const std::size_t possible = static_cast<std::size_t>(1 + dam_count +
      sire_count) + static_cast<std::size_t>(dam_count) * sire_count;
    configurations.reserve(possible);
    auto add_configuration = [&](int dam, int sire) {

      ConfigurationScore configuration;
      configuration.dam = dam;
      configuration.sire = sire;
      configuration.nparents = (dam > 0 ? 1 : 0) + (sire > 0 ? 1 : 0);
      ConfigurationScore incumbent;
      incumbent.dam = starting_dam; incumbent.sire = starting_sire;
      incumbent.nparents = (starting_dam > 0) + (starting_sire > 0);
      set_incumbent_policy(input, incumbent, configuration);

      const auto proposal_started = input.bounded_kernel_profile ?
        std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
      const int proposal_reason = state.proposal_reason(child, dam, sire);

      configuration.hard_feasible = proposal_reason == 0;
      if (!configuration.hard_feasible) {
        configuration.status = 1;
        configuration.reason = proposal_reason;
        configurations.push_back(configuration);
        return;
      }
      if (configuration.nparents == 2) {
        const auto trio_started = (input.benchmark_opt5_diag_mode || input.bounded_kernel_profile) ?
          std::chrono::steady_clock::now() :
          std::chrono::steady_clock::time_point();
        configuration.trio_me = count_trio_me(
          compiled, child, dam, sire);

        configuration.trio_me_pass = configuration.trio_me < 0 ||
          configuration.trio_me <= input.max_mendelian_errors;
        if (configuration.trio_me_pass && input.trio_lr_filter_ln > R_NegInf &&
            !(dam == starting_dam && sire == starting_sire))
          configuration.trio_me_pass =
            trio_llr_ln(compiled, child, dam, sire) >= input.trio_lr_filter_ln;
      }
      configurations.push_back(configuration);
    };
    // Relative duo screen: a candidate whose unconditional parent-offspring
    // likelihood ratio trails its side's best by more than the margin is
    // not proposed, alone or in a pair. Current parents are always proposed.
    std::vector<unsigned char> dam_screened(static_cast<std::size_t>(dam_count), 0U),
      sire_screened(static_cast<std::size_t>(sire_count), 0U);
    if (input.duo_lr_margin_ln > 0.0 && (!input.duo_lr_margin_recheck_only || input.recheck_only)) {
      auto screen = [&](int begin, int end, const std::vector<int>& ids, int current, bool dam_side,
                        std::vector<unsigned char>& out) {
        std::vector<double> lr(static_cast<std::size_t>(end - begin));
        double best = R_NegInf;
        for (int i = begin; i < end; ++i) {
          lr[i - begin] = dam_side ? trio_llr_ln(compiled, child, ids[i], 0) : trio_llr_ln(compiled, child, 0, ids[i]);
          best = std::max(best, lr[i - begin]);
        }
        for (int i = begin; i < end; ++i)
          out[i - begin] = ids[i] != current && lr[i - begin] < best - input.duo_lr_margin_ln;
      };
      screen(dam_begin, dam_end, input.dam_ids, starting_dam, true, dam_screened);
      screen(sire_begin, sire_end, input.sire_ids, starting_sire, false, sire_screened);
    }
    add_configuration(0, 0);
    for (int i = dam_begin; i < dam_end; ++i)
      if (!dam_screened[i - dam_begin]) add_configuration(input.dam_ids[i], 0);
    for (int i = sire_begin; i < sire_end; ++i)
      if (!sire_screened[i - sire_begin]) add_configuration(0, input.sire_ids[i]);
    const bool pair_rescue_active=input.pair_rescue_enabled||input.pair_rescue_shadow;
    pap_internal::pair_rescue::Visit pair_gate;
    pair_gate.round=input.progress_parentage_round;pair_gate.focal=child;
    pair_gate.revision=graph_before;pair_gate.recheck=input.recheck_only;
    pair_gate.completion_round=identity_completion_round;
    pair_gate.completion_rows=int(completed_identity_competitors.size());
    std::vector<pap_internal::pair_rescue::Candidate> pair_dams,pair_sires;
    std::vector<int> pair_dam_index,pair_sire_index;
    if(pair_rescue_active){
      const auto started=std::chrono::steady_clock::now();
      pair_gate.configurations_before=int(configurations.size());
      pair_dam_index.assign(input.n+1,-1);pair_sire_index.assign(input.n+1,-1);
      auto* live=pap_internal::native_demand::live;
      for(const auto& row:configurations)if(row.nparents==1){
        const int side=row.dam?1:2,id=row.dam?row.dam:row.sire;
        pap_internal::pair_rescue::Candidate c;c.id=id;
        c.feasible=row.hard_feasible&&row.status==0;
        c.incumbent=id==(side==1?starting_dam:starting_sire);
        if(c.feasible){c.score=live->single_candidate(id,side,input.Tassign_ln-input.pair_rescue_margin_ln,input.pair_rescue_shadow,input.pair_rescue_support_only);++pair_gate.singles;
          pair_gate.support_shortcuts+=c.score.relationship_skipped&&c.score.support_rejected;
          pair_gate.cache_hits+=c.score.cached;pair_gate.unknown+=pap_internal::pair_rescue::unknown(c.score);}
        auto& list=side==1?pair_dams:pair_sires;
        (side==1?pair_dam_index:pair_sire_index)[id]=int(list.size());list.push_back(c);
      }
      pap_internal::pair_rescue::select(pair_dams,input.Tassign_ln,input.pair_rescue_margin_ln,input.pair_rescue_top_k);
      pap_internal::pair_rescue::select(pair_sires,input.Tassign_ln,input.pair_rescue_margin_ln,input.pair_rescue_top_k);
      for(int side:{1,2})for(const auto& c:side==1?pair_dams:pair_sires){
        pair_gate.passing+=c.passing;pair_gate.retained+=c.retained;
        if(c.feasible&&!c.retained)++(side==1?pair_gate.rank_capped_dams:pair_gate.rank_capped_sires);
        if(c.passing&&!pap_internal::pair_rescue::unknown(c.score))++(side==1?pair_gate.dam_floor_pass:pair_gate.sire_floor_pass);
        if(input.pair_rescue_diagnostics)compiled.pair_rescue_scores.push_back(
          {pair_gate.round,child,graph_before,int(input.recheck_only),side,c});
      }
      if(input.pair_rescue_enabled&&input.pair_rescue_support_only&&input.config->likelihood_ancestral){
        configurations.erase(std::remove_if(configurations.begin(),configurations.end(),[&](const auto& row){
          if(row.nparents!=1)return false;
          const auto& c=row.dam?pair_dams.at(pair_dam_index[row.dam]):pair_sires.at(pair_sire_index[row.sire]);
          return !c.retained&&!completed_identity_competitors.count({row.dam,row.sire});
        }),configurations.end());
      }
      pair_gate.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    }
    for (int i = dam_begin; i < dam_end; ++i)
      for (int j = sire_begin; j < sire_end; ++j)
        if (input.dam_ids[i] != input.sire_ids[j] && !dam_screened[i - dam_begin] && !sire_screened[j - sire_begin]) {
          const int dam=input.dam_ids[i],sire=input.sire_ids[j];
          bool candidate_keep=true;
          const bool current_pair=dam==starting_dam&&sire==starting_sire;
          if(pair_rescue_active){
            ++pair_gate.pairs;
            const auto& d=pair_dams.at(pair_dam_index[dam]);const auto& s=pair_sires.at(pair_sire_index[sire]);
            candidate_keep=pap_internal::pair_rescue::pair(d,s,current_pair)||
              completed_identity_competitors.count({dam,sire});
            if(!candidate_keep){if(!d.passing&&!s.passing)++pair_gate.failed_both;else ++pair_gate.rank_pruned;}
            if(input.pair_rescue_enabled&&!candidate_keep)continue;
          }
          add_configuration(dam,sire);
          if(pair_rescue_active&&(candidate_keep||input.pair_rescue_diagnostics)){
            const auto& row=configurations.back();
            // Structural eligibility precedes even the cheap trio likelihood.
            const double trio_score=row.hard_feasible&&row.trio_me_pass?
              trio_llr_ln(compiled,child,dam,sire):NA_REAL;
            const bool trio_keep=current_pair||(!input.completion_rank_only&&completed_identity_competitors.count({dam,sire}))||!row.hard_feasible||
              (row.trio_me_pass&&(std::isnan(trio_score)||trio_score>=input.pair_rescue_trio_ln));
            if(input.pair_rescue_diagnostics)compiled.pair_rescue_pairs.push_back(
              {pair_gate.round,child,graph_before,int(input.recheck_only),dam,sire,current_pair,
               row.hard_feasible,row.trio_me_pass,candidate_keep,trio_score});
            if(candidate_keep){
              if(trio_keep)++pair_gate.kept;
              else{++pair_gate.trio_pruned;if(input.pair_rescue_enabled)configurations.pop_back();}
            }
          }
        }
    if(pair_rescue_active){
      pair_gate.configurations_before+=pair_gate.pairs;pair_gate.configurations_after=int(configurations.size());
      if(input.pair_rescue_enabled&&input.pair_rescue_support_only&&input.config->likelihood_ancestral){
        std::vector<int> roots;
        for(const auto& row:configurations)if(row.hard_feasible&&row.status==0)
          for(int id:{row.dam,row.sire})if(id)roots.push_back(id);
        // Completion retries must not change a frozen union after numerical
        // scoring. Bind every possible completion endpoint up front. The
        // default configuration-local model still reads only its own core.
        if(input.complete_identity_competitors){
          roots.insert(roots.end(),input.dam_ids.begin()+dam_begin,input.dam_ids.begin()+dam_end);
          roots.insert(roots.end(),input.sire_ids.begin()+sire_begin,input.sire_ids.begin()+sire_end);
        }
        // Keep incumbent observations in the same fixed retention comparison.
        for(int id:{starting_dam,starting_sire})if(id)roots.push_back(id);
        auto counts=pap_internal::native_demand::live->freeze_identity_roots(roots);
        pair_gate.roots_before=counts[0];pair_gate.roots_after=counts[1];
      }
      compiled.pair_rescue_visits.push_back(pair_gate);
      pap_internal::pair_rescue::completion_work=&compiled.pair_rescue_visits.back();
    }
    // Relative trio screen: a candidate pair whose unconditional trio
    // likelihood ratio trails the visit's best pair by more than the margin
    // is screened like a Mendelian-error failure. The incumbent pair stays.
    if (input.trio_lr_margin_ln > 0.0 && (!input.trio_lr_margin_recheck_only || input.recheck_only)) {
      std::vector<double> pair_lr(configurations.size(), R_NegInf);
      double best_pair_lr = R_NegInf;
      for (std::size_t index = 0; index < configurations.size(); ++index) {
        const auto& row = configurations[index];
        if (row.nparents != 2 || !row.hard_feasible || !row.trio_me_pass) continue;
        pair_lr[index] = trio_llr_ln(compiled, child, row.dam, row.sire);
        best_pair_lr = std::max(best_pair_lr, pair_lr[index]);
      }
      for (std::size_t index = 0; index < configurations.size(); ++index) {
        auto& row = configurations[index];
        if (row.nparents != 2 || !row.hard_feasible || !row.trio_me_pass) continue;
        if (row.dam == starting_dam && row.sire == starting_sire) continue;
        if (pair_lr[index] < best_pair_lr - input.trio_lr_margin_ln) { row.trio_me_pass = false;  }
      }
    }

    // The external storage is created after the fixed candidate configuration
    // vector is complete so its element addresses are stable throughout the
    // locus-major pass.  Both modes use its ordinary members; only full adds
    // compound/contextual members.
    std::vector<SideTotals> full_relationship_members;
    if (!input.full_relationship_delta.empty()) {
      full_relationship_members.resize(configurations.size());
      for (std::size_t index = 0; index < configurations.size(); ++index)
        configurations[index].full_relationship_members =
          &full_relationship_members[index];
    }

    int single_configurations = 0;
    int pair_configurations = 0;
    for (const ConfigurationScore& configuration : configurations) {
      if (configuration.nparents == 1) ++single_configurations;
      if (configuration.nparents == 2) ++pair_configurations;
    }

    const auto message_trial_failed = [&]() {
      return input.message_scoring && !input.private_complete_message_trial &&
        progress_telemetry.message_partition_unavailable > message_declines_before;
    };
    const auto abandon_message_trial = [&]() {

        message_fallback = true;

    };
    const bool capture_configuration_loci =
      input.configuration_locus_diagnostics && diagnostic_focal &&
      (input.configuration_locus_focal_position == -1 ||
       static_cast<int>(focal_position + 1) ==
         input.configuration_locus_focal_position);
    std::vector<double> captured_locus_scores;
    if (capture_configuration_loci)
      captured_locus_scores.assign(
        configurations.size() * static_cast<std::size_t>(compiled.loci), NA_REAL);
    // Freeze the raw scientific request from the complete candidate universe
    // before any score-dependent work elimination or relationship construction.



    std::vector<RelationshipContexts> side_contexts(configurations.size());
    if (state.pass.bounded_messages && !input.elide_raw()) {
      const auto raw_boundary = build_raw_common_boundary(compiled, state, child, configurations);
      if(auto* live=pap_internal::native_demand::live;live && raw_boundary)live->observe_raw_catalog(child,*raw_boundary);
      for (std::size_t index = 0; index < side_contexts.size(); ++index) {
        side_contexts[index].raw_boundary = raw_boundary;
        side_contexts[index].raw_request_index = index;
      }
    }
    // Reference and exhaustive observers retain the complete old ordering.
    const bool stage_raw_support = input.lazy_independent_gates &&
      state.pass.bounded_messages && compiled.loci > 0 &&
      !input.configuration_diagnostics && !input.configuration_locus_diagnostics &&
      !input.benchmark_baseline_mode &&
      input.config->evaluation.relationship_common_scaled_fast &&
      !input.private_test_relationship_comparison_exception &&
      !input.private_test_raw_comparison_exception &&
      !input.private_test_configuration_post_score_exception;
    worker_team.configure_raw_arithmetic(stage_raw_support && input.raw_scaled_fast);
    state.pass.raw_scaled_fast = stage_raw_support && input.raw_scaled_fast;
    double raw_support_scoring_seconds = 0.0;
    if (stage_raw_support) {

      worker_team.probe_required_raw_backend(input, side_contexts, configurations);
      if (message_trial_failed()) {

        abandon_message_trial(); goto retry_frozen_focal_model;
      }
      const auto raw_started = input.phase_timing ?
        std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();

      worker_team.score_raw_panel(input, child, side_contexts, configurations);
      if (message_trial_failed()) { abandon_message_trial(); goto retry_frozen_focal_model; }
      if (input.phase_timing)
        raw_support_scoring_seconds = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - raw_started).count();
    }
    if(stage_raw_support&&!window_widened&&
       window_requires_unknown_identity(input,configurations)){
      window_widened=true;input.gate_window_ln=0;

      // The pending memo key describes the original window. As for a late
      // widened retry, do not publish a result under that execution key.
      mark_raw_support_rejections(input,compiled.loci,configurations);
    }
    // Derive retained contexts only for configurations whose independent gates
    // still need a relationship answer. Positive full rows retain every direct
    // and compound parental numerator and all applicable ordinary alternatives.
    if (!input.full_relationship_delta.empty())
      for (std::size_t configuration_index = 0;
           configuration_index < configurations.size(); ++configuration_index) {
        const ConfigurationScore& configuration = configurations[configuration_index];
        if (!configuration.hard_feasible || configuration.status != 0 ||
            !relationship_stage_needed(configuration)) continue;
        if (input.lazy_independent_gates && !input.configuration_diagnostics &&
            !input.configuration_locus_diagnostics &&
            (!configuration.trio_me_pass)) continue;
        RelationshipContexts prepared = build_side_contexts(
          state, child, configuration, configurations);
        prepared.raw_boundary = side_contexts[configuration_index].raw_boundary;
        prepared.raw_request_index = side_contexts[configuration_index].raw_request_index;
        side_contexts[configuration_index] = std::move(prepared);
      }
    std::vector<unsigned char> numerical_cached(configurations.size(), 0U);



    // The bounded locus interval is deliberately only a cheap C++ check.  A
    // C++ -> R callback is wall-clock throttled so the benchmark baseline does
    // not turn every five completed loci into file I/O and allocation work.
    auto last_locus_heartbeat = std::chrono::steady_clock::now();

    LocusResult& locus_result = worker.locus_results.front();
    const bool compact_native_results=pap_internal::native_demand::live&&pap_internal::native_demand::live->engine_compact_results&&!input.benchmark_baseline_mode&&!input.message_scoring&&!capture_configuration_loci;
    if (input.benchmark_baseline_mode)
      locus_result.reset(worker.owned_loci.front(), configurations.size());

    // Score locus-major while preserving the per-configuration locus order.
    // Each objective still accumulates locus 1, 2, ... exactly as before, but
    // cavity/factor memo state for a completed locus can be released instead
    // of retaining equivalent recursion forests for every locus until the
    // focal decision is complete.
    for (std::size_t index = 0U; index < configurations.size(); ++index)
      if (configurations[index].hard_feasible && !configurations[index].raw_precomputed &&
          !numerical_cached[index]) configurations[index].objective_ln = 0.0;
    const auto locus_scoring_started = input.phase_timing ?
      std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
    if (input.benchmark_baseline_mode) {
      for (std::size_t configuration_index = 0;
           configuration_index < configurations.size(); ++configuration_index) {
        ConfigurationScore& configuration = configurations[configuration_index];
        if (!configuration.hard_feasible || configuration.status != 0) continue;
        for (std::size_t locus_index = 0;
             locus_index < workspaces.size(); ++locus_index) {
          if (configuration.status != 0) break;
          locus_result.locus_index = static_cast<int>(locus_index);
          ConfigurationResult& result =
            locus_result.configurations[configuration_index];
          evaluate_configuration_locus(
            input, worker.pass, *workspaces[locus_index], locus_index, child,
            configuration,
            side_contexts[configuration_index], result);
          if (capture_configuration_loci)
            captured_locus_scores[configuration_index *
                                  static_cast<std::size_t>(compiled.loci) +
                                  locus_index] = result.locus_objective_ln;
          reduce_configuration_locus(
            result, configuration);
        }
      }

    } else {
      for (std::size_t locus_index = 0;
           locus_index < workspaces.size(); ++locus_index) {
        LocusWorkspace& workspace = *workspaces[locus_index];
        locus_result.reset(static_cast<int>(locus_index), compact_native_results?1:configurations.size());
        if(compact_native_results)pap_internal::native_demand::live->compact_packets_avoided+=configurations.size()-1;
        for (std::size_t configuration_index = 0;
             configuration_index < configurations.size();
             ++configuration_index) {
          const ConfigurationScore& configuration =
            configurations[configuration_index];
          if (numerical_cached[configuration_index] || !configuration.hard_feasible ||
              configuration.status != 0 || !relationship_stage_needed(configuration)) continue;

            evaluate_configuration_locus(
              input, worker.pass, workspace, locus_index, child, configuration,
              side_contexts[configuration_index],
              locus_result.configurations[compact_native_results?0:configuration_index]);

          if (message_trial_failed()) break;
          if(compact_native_results)reduce_configuration_locus(locus_result.configurations[0],configurations[configuration_index]);
        }
        if (message_trial_failed()) { abandon_message_trial(); goto retry_frozen_focal_model; }
        for (std::size_t configuration_index = 0;
             !compact_native_results && configuration_index < configurations.size();
             ++configuration_index) {
          ConfigurationResult& result =
            locus_result.configurations[compact_native_results?0:configuration_index];
          if (!result.scored) continue;
          ConfigurationScore& configuration = configurations[configuration_index];
          if (capture_configuration_loci)
            captured_locus_scores[configuration_index *
                                  static_cast<std::size_t>(compiled.loci) +
                                  locus_index] = result.locus_objective_ln;
          reduce_configuration_locus(
            result, configuration);

        }

        const bool final_locus = locus_index + 1 == workspaces.size();
        const bool checkpoint_due =
          (static_cast<int>(locus_index) + 1) %
            input.progress_locus_interval == 0 || final_locus;
        if (input.progress_locus_heartbeats && checkpoint_due && (final_locus ||
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
              last_locus_heartbeat).count() >= input.progress_heartbeat_seconds)) {

          last_locus_heartbeat = std::chrono::steady_clock::now();
        }


      }
    }
    // A complete cached panel may replace prepass Q. Its numerical bound
    // follows that value; redo numeric protection against the current universe.
    worker_team.refine_raw_panels(input, side_contexts, configurations);
    for (auto& configuration : configurations)
      restore_cached_raw_parent_alias(input, static_cast<std::size_t>(compiled.loci),
        configuration, configuration.full_relationship_members);

    worker_team.refine_relationship_panels(input, child, side_contexts, configurations);
    if (message_trial_failed()) { abandon_message_trial(); goto retry_frozen_focal_model; }
    worker_team.capture_result_storage(input);

    if (input.message_scoring &&
        progress_telemetry.message_partition_unavailable > message_declines_before) {
      // Nothing from this trial reaches a chooser, mutation or returned
      // comparison. Repeat the same frozen focal under the retained model,
      // including all raw identities and applicable ordinary/full rivals.
      // Cache keys include the scorer model; never relabel mixed operands.

      message_fallback = true;

      goto retry_frozen_focal_model;
    }

     // uncached preparation and scoring; all decisions below are current.

    const auto decision_started = input.phase_timing ?
      std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
    for (const ConfigurationScore& configuration : configurations) {
      if (!configuration.hard_feasible) continue;
      if (configuration.status == 2)
        {}
      else {}
    }

    auto find_configuration = [&](int dam, int sire) {
      for (std::size_t i = 0; i < configurations.size(); ++i)
        if (configurations[i].dam == dam && configurations[i].sire == sire)
          return static_cast<int>(i);
      return -1;
    };
    const int current_index = find_configuration(starting_dam, starting_sire);
    if (current_index < 0)
      engine_stop("whole-pass candidate arrays omitted the current parent configuration");
    const auto diagnostics_started = input.bounded_kernel_profile ?
      std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
    compute_diagnostics(input, current_index, configurations,
      nullptr,
      nullptr);

    int single_total = 0, pair_total = 0;
    int single_hard = 0, pair_hard = 0;
    int single_me = 0, pair_me = 0;
    int single_relationship = 0, pair_relationship = 0;
    int single_eligible = 0, pair_eligible = 0;
    for (const ConfigurationScore& configuration : configurations) {
      if (configuration.nparents == 1) {
        ++single_total;
        if (configuration.hard_feasible) ++single_hard;
        if (configuration.hard_feasible && configuration.trio_me_pass) ++single_me;
        if (configuration.hard_feasible &&
            configuration.relationship_pass) ++single_relationship;
        if (configuration.eligible) ++single_eligible;
      } else if (configuration.nparents == 2) {
        ++pair_total;
        if (configuration.hard_feasible) ++pair_hard;
        if (configuration.hard_feasible && configuration.trio_me_pass) ++pair_me;
        if (configuration.hard_feasible && configuration.trio_me_pass &&
            configuration.relationship_pass)
          ++pair_relationship;
        if (configuration.eligible) ++pair_eligible;
      }
    }
    bool near_tie = false;
    int rejection_reason = 1;
    const auto choose_started =
      (input.phase_timing || input.bounded_kernel_profile) ?
      std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
    const int chosen_index = choose_configuration(
      input, current_index, configurations, near_tie, rejection_reason);
    if (input.gate_window_ln > 0.0 && !window_widened) {
      bool eligible = false, deferred = false;
      int best = -1;
      for (std::size_t i = 0; i < configurations.size(); ++i) {
        const auto& configuration = configurations[i];
        if (configuration.raw_comparison_available && R_finite(configuration.objective_ln) &&
            (best < 0 || configuration.objective_ln > configurations[best].objective_ln)) best = static_cast<int>(i);
        if (configuration.nparents == 0) continue;
        eligible = eligible || configuration.eligible;
        deferred = deferred || (configuration.gate_deferred &&
          (!input.joint_identity_objective||identity_competitor(configuration)));
      }
      if (!eligible && deferred && best >= 0 && configurations[best].nparents > 0) {
        window_widened = true;
        input.gate_window_ln = 0.0;

        goto retry_frozen_focal_model;
      }
    }
    if(input.complete_identity_competitors&&input.pair_rescue_enabled&&
       input.config->likelihood_ancestral&&input.joint_identity_objective&&chosen_index>=0&&
       (!input.completion_rank_only||pair_gate.rank_pruned>0)){
      const auto& provisional=configurations[chosen_index];
      if(provisional.nparents>0){
        bool added=false;
        const auto peers=pap_internal::pair_rescue::completions(
          provisional.dam,provisional.sire,
          std::span<const int>(input.dam_ids).subspan(dam_begin,dam_count),
          std::span<const int>(input.sire_ids).subspan(sire_begin,sire_count));
        for(const auto& peer:peers)if(find_configuration(peer.first,peer.second)<0){
          if(input.completion_rank_only){
            if(!peer.first||!peer.second)continue;
            const int di=pair_dam_index[peer.first],si=pair_sire_index[peer.second];
            if(di<0||si<0||!pap_internal::pair_rescue::rank_omitted_pair(pair_dams[di],pair_sires[si]))continue;
            // The rank cap alone is relaxed. Preserve structural, Mendelian
            // and cheap trio screens before adding an expensive context.
            if(state.proposal_reason(child,peer.first,peer.second)!=0)continue;
            if(count_trio_me(compiled,child,peer.first,peer.second)>input.max_mendelian_errors)continue;
            const double trio_score=trio_llr_ln(compiled,child,peer.first,peer.second);
            if(!std::isnan(trio_score)&&trio_score<input.pair_rescue_trio_ln)continue;
          }
          added|=completed_identity_competitors.insert(peer).second;
        }
        // Each retry admits at least one new finite candidate identity. The
        // graph, thresholds, kernels and cut contexts remain frozen. Unknown
        // comparisons retain the ordinary chooser's uncertainty semantics.
        if(added){++identity_completion_round;goto retry_frozen_focal_model;}
      }
    }
    // The chooser has finished, and the graph has not moved. Probe extra
    // identities without replacing any scored row or repeating the chooser.
    if(input.raw_joint_census&&input.config->likelihood_ancestral){
      auto* live=pap_internal::native_demand::live;
      if(!live)engine_stop("raw/joint census requires the M_ped2 live scorer");
      const auto started=std::chrono::steady_clock::now();
      std::array<double,14> r{};r[0]=pair_gate.round;r[1]=child;r[2]=graph_before;r[3]=input.recheck_only;
      double null=NA_REAL;
      for(const auto& c:configurations)if(!c.nparents&&c.raw_comparison_available)null=c.objective_ln;
      r[12]=null;double low=R_finite(null)?null:INFINITY,high=R_finite(null)?null:-INFINITY;
      for(const auto& c:configurations){
        if(!c.hard_feasible||c.status!=0||!c.trio_me_pass)continue;
        ++r[4];r[5]+=c.nparents==2;
        if(!c.raw_comparison_available||c.raw_loci_accumulated!=size_t(compiled.loci)||!R_finite(c.objective_ln)){++r[6];continue;}
        const double joint=c.nparents?live->joint_identity(c.dam,c.sire):0.;
        if(!R_finite(joint)){++r[7];r[8]=INFINITY;r[10]=c.dam;r[11]=c.sire;continue;}
        const double residual=c.objective_ln-joint;
        low=std::min(low,residual);high=std::max(high,residual);
        const double error=R_finite(null)?std::abs(residual-null):NA_REAL;
        if(error>r[8]){r[8]=error;r[10]=c.dam;r[11]=c.sire;}
      }
      r[9]=R_finite(low)&&R_finite(high)?high-low:NA_REAL;
      r[13]=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
      compiled.raw_joint_checks.push_back(r);
    }
    if(input.pair_rescue_calibration){
      namespace rescue=pap_internal::pair_rescue;
      for(std::size_t index=0;index<configurations.size();++index){
        const auto& row=configurations[index];
        if(row.nparents!=2||!identity_competitor(row))continue;
        const auto& d=pair_dams.at(pair_dam_index[row.dam]);
        const auto& s=pair_sires.at(pair_sire_index[row.sire]);
        const bool current=row.dam==starting_dam&&row.sire==starting_sire;
        const double cheap=trio_llr_ln(compiled,child,row.dam,row.sire);
        compiled.pair_rescue_requirements.push_back({pair_gate.round,child,graph_before,
          int(input.recheck_only),row.dam,row.sire,rescue::required_rank(d,pair_dams),
          rescue::required_rank(s,pair_sires),rescue::required_margin(d,s,input.Tassign_ln),
          current,int(index)==chosen_index,current||std::isnan(cheap)||cheap>=input.pair_rescue_trio_ln});
      }
    }
    for (const auto& configuration : configurations)

    {}

    const double chooser_seconds = input.phase_timing ? std::chrono::duration<double>(std::chrono::steady_clock::now()-choose_started).count() : 0.;
    ConfigurationScore& current = configurations[current_index];
    ConfigurationScore& selected = configurations[chosen_index];
    if(diagnostic_sequence)for(std::size_t config=0;config<configurations.size();++config) {
      const auto& c=configurations[config];DiagnosticPanelRow r;
      r.visit=visit_index+1;r.focal=child;r.graph=graph_before;r.dam=c.dam;r.sire=c.sire;
      r.objective=c.objective_ln;r.dam_support=c.dam_support_ln;r.sire_support=c.sire_support_ln;
      r.pair_support=c.pair_support_ln;r.full_relationship=c.full_relationship_ln;r.identity=c.identity_margin_ln;
      r.joint=decision_objective(c);r.competitor=identity_competitor(c);
      r.raw_available=c.raw_comparison_available;r.full_available=c.relationship_comparison_available;
      if(c.full_relationship_members){const auto& members=*c.full_relationship_members;auto reduced=full_reduce_members(&members,c.nparents,NA_REAL);
       if(reduced.comparison_available&&reduced.selected_side>=0&&members.relationship_context[reduced.selected_side].common_boundary){const auto& boundary=*members.relationship_context[reduced.selected_side].common_boundary;
        r.strongest_rival="side"+std::to_string(reduced.selected_side+1)+":";
        for(int id:reduced.alternative_ties){if(r.strongest_rival.back()!=':')r.strongest_rival+=';';
         r.strongest_rival+=id>=0&&id<25?boundary.models[id].topology.name:id>=25&&size_t(id-25)<boundary.full_boundary->members.size()?boundary.full_boundary->members[id-25].name:"unknown";}
       }}
      r.current=config==static_cast<std::size_t>(current_index);r.selected=config==static_cast<std::size_t>(chosen_index);r.eligible=c.eligible;
      diagnostic_panels.push_back(r);
    }

    if (input.configuration_diagnostics && diagnostic_focal) {
      configuration_diagnostics.reserve(
        configuration_diagnostics.size() + configurations.size());
      for (std::size_t configuration_position = 0;
           configuration_position < configurations.size();
           ++configuration_position) {
        const ConfigurationScore& configuration =
          configurations[configuration_position];
        ConfigurationDiagnostic record;
        record.focal_position = static_cast<int>(focal_position + 1);
        record.focal = child;
        record.graph_revision = graph_before;
        record.current_dam_before_test = state.pass.dam[child];
        record.current_sire_before_test = state.pass.sire[child];
        record.configuration_position = static_cast<int>(
          configuration_position + 1);
        record.dam = configuration.dam;
        record.sire = configuration.sire;
        record.nparents = configuration.nparents;
        record.status = configuration.status;
        record.reason = configuration.reason;
        record.route = configuration.route;
        record.trio_me = configuration.trio_me;
        record.objective_ln = configuration.objective_ln;
        record.identity_objective_ln = decision_objective(configuration);
        record.raw_reference_error = configuration.raw_reference_error;
        record.raw_comparison_available = configuration.raw_comparison_available ? 1 : 0;
        record.raw_status = configuration.raw_status;
        record.dam_support_ln = configuration.dam_support_ln;
        record.sire_support_ln = configuration.sire_support_ln;
        record.pair_support_ln = configuration.pair_support_ln;
        record.identity_margin_ln = configuration.identity_margin_ln;
        record.hard_feasible = configuration.hard_feasible ? 1 : 0;
        record.hard_valid = configuration.hard_valid ? 1 : 0;
        record.trio_me_pass = configuration.trio_me_pass ? 1 : 0;
        record.relationship_pass = configuration.relationship_pass ? 1 : 0;
        record.relationship_comparison_available =
          configuration.relationship_comparison_available ? 1 : 0;
        record.full_relationship_gate_pass =
          configuration.full_relationship_gate_pass ? 1 : 0;
        record.full_relationship_ln = configuration.full_relationship_ln;
        record.full_operational_ln = configuration.full_operational_ln;
        record.full_alternative_ln = configuration.full_alternative_ln;
        record.assignment_confident = configuration.assignment_confident ? 1 : 0;
        record.eligible_for_retention =
          configuration.eligible_for_retention ? 1 : 0;
        record.eligible_for_addition =
          configuration.eligible_for_addition ? 1 : 0;
        record.eligible_for_replacement =
          configuration.eligible_for_replacement ? 1 : 0;
        record.eligible_for_removal =
          configuration.eligible_for_removal ? 1 : 0;
        record.eligible = configuration.eligible ? 1 : 0;
        record.near_tie = configuration.near_tie ? 1 : 0;
        record.is_current = static_cast<int>(configuration_position) ==
          current_index ? 1 : 0;
        record.is_selected = static_cast<int>(configuration_position) ==
          chosen_index ? 1 : 0;
        record.focal_rejection_reason = rejection_reason;
        if (configuration.full_relationship_members != nullptr)
          record.full_relationship_member_ln =
            configuration.full_relationship_members->value;
        if (configuration.full_relationship_members != nullptr)
          record.full_parent_numerator_ln =
            configuration.full_relationship_members->parent_value;
        if (configuration.full_relationship_members != nullptr)
          record.full_contextual_member_ln =
            configuration.full_relationship_members->contextual_value;
        if (configuration.full_relationship_members != nullptr) {
          const SideTotals& members =
            *configuration.full_relationship_members;
          record.side_reference_uncertainty = members.relationship_common_error_bound;
          if (members.full_common_requested[0] || members.full_common_requested[1]) {
            record.full_common_boundary = 1;
            const SideReduction reduced = full_reduce_members(
              &members, configuration.nparents, NA_REAL);
            record.full_limiting_side = reduced.selected_side;
            record.full_side_relationship_ln = reduced.side_relationship_ln;
            for (std::size_t side = 0U; side < 2U; ++side)
              record.side_reference_uncertainty[side] = std::max(
                record.side_reference_uncertainty[side], reduced.side_reference_uncertainty[side]);
            if (reduced.comparison_available && reduced.selected_side >= 0) {
              const auto& boundary = *members.relationship_context[
                static_cast<std::size_t>(reduced.selected_side)].common_boundary;
              const auto names = [&](const std::vector<int>& ids) {
                std::string joined;
                for (int id : ids) {
                  const std::string name = id >= 0 && id < 25 ?
                    boundary.models[static_cast<std::size_t>(id)].topology.name :
                    id >= 25 && static_cast<std::size_t>(id - 25) <
                      boundary.full_boundary->members.size() ?
                    boundary.full_boundary->members[static_cast<std::size_t>(id - 25)].name :
                    "unknown_model";
                  if (!joined.empty()) joined += ";";
                  joined += name;
                }
                return joined;
              };
              record.full_parent_model = names(reduced.operational_ties);
              record.full_alternative_model = names(reduced.alternative_ties);
            }
          }
          for (int side = 0; side < 2; ++side) {
            if (!members.relationship_context_seen[static_cast<std::size_t>(side)])
              continue;
            const RelationshipContext& context =
              members.relationship_context[static_cast<std::size_t>(side)];
            RelationshipProvenanceDiagnostic& provenance =
              record.relationship_provenance[static_cast<std::size_t>(side)];
            provenance.focal = context.focal;
            provenance.candidate = context.candidate;
            provenance.tested_parent_side = context.proposed_parent_side;
            provenance.proposed_dam = context.proposed_dam;
            provenance.proposed_sire = context.proposed_sire;
            provenance.focal_dam = context.focal_parents.dam;
            provenance.focal_sire = context.focal_parents.sire;
            provenance.candidate_dam = context.candidate_parents.dam;
            provenance.candidate_sire = context.candidate_parents.sire;
            provenance.fixed_other_parent = context.proposed_parent_side == 1 ?
              context.proposed_sire : (context.proposed_parent_side == 2 ?
              context.proposed_dam : 0);
            provenance.has_relevant_context = context.has_relevant_context ? 1 : 0;
            provenance.focal_endpoint_excluded_child =
              context.focal_endpoint_excluded_child;
            provenance.candidate_endpoint_excluded_child =
              context.candidate_endpoint_excluded_child;
            provenance.fs_orientation_mask = context.fs_orientation_mask;
            provenance.hs_orientation_mask = context.hs_orientation_mask;
            provenance.gp_orientation_mask = context.gp_orientation_mask;
            provenance.fa_orientation_mask = context.fa_orientation_mask;
            provenance.ha_orientation_mask = context.ha_orientation_mask;
            provenance.loci_accumulated = members.relationship_loci_accumulated[
              static_cast<std::size_t>(side)];
            provenance.context_signature = context.signature;
            const std::size_t member_base = static_cast<std::size_t>(
              side * kRelationshipCount);
            for (int relationship = 0;
                 relationship < kRelationshipCount; ++relationship)
              provenance.partition_ln[static_cast<std::size_t>(relationship)] =
                members.relationship_partition_value[member_base +
                  static_cast<std::size_t>(relationship)];
            provenance.proposed_po_partition_ln =
              members.relationship_po_partition_value[static_cast<std::size_t>(side)];
            for (int orientation = 0; orientation < 22; ++orientation)
              provenance.orientation_partition_ln[static_cast<std::size_t>(orientation)] =
                members.relationship_orientation_partition_value[
                  static_cast<std::size_t>(side * 22 + orientation)];
            record.relationship_provenance_seen[static_cast<std::size_t>(side)] = 1;
          }
        }

        configuration_diagnostics.push_back(record);
      }
    }
    if (capture_configuration_loci) {
      for (std::size_t configuration_index = 0;
           configuration_index < configurations.size(); ++configuration_index)
        for (std::size_t locus_index = 0;
             locus_index < static_cast<std::size_t>(compiled.loci); ++locus_index) {
          ConfigurationLocusDiagnostic record;
          record.focal_position = static_cast<int>(focal_position + 1);
          record.focal = child;
          record.graph_revision = graph_before;
          record.configuration_position = static_cast<int>(configuration_index + 1);
          record.locus = static_cast<int>(locus_index + 1);
          record.objective_ln = captured_locus_scores[
            configuration_index * static_cast<std::size_t>(compiled.loci) +
            locus_index];
          configuration_locus_diagnostics.push_back(record);
        }
    }

    // Observe only after chooser/gate evaluation is complete. The diagnostic
    // cannot supply likelihoods to any decision made during this visit.

    VisitSummary summary;summary.chooser_seconds=chooser_seconds;
    summary.focal = child;
    summary.starting_dam = starting_dam;
    summary.starting_sire = starting_sire;
    summary.selected_dam = selected.dam;
    summary.selected_sire = selected.sire;
    summary.current_objective_ln = input.joint_identity_objective ? decision_objective(current) : current.objective_ln;
    summary.selected_objective_ln = input.joint_identity_objective ? decision_objective(selected) : selected.objective_ln;
    summary.dam_support_ln = selected.dam_support_ln;
    summary.sire_support_ln = selected.sire_support_ln;
    summary.pair_support_ln = selected.pair_support_ln;
    summary.identity_margin_ln = selected.identity_margin_ln;
    summary.full_relationship_ln = selected.full_relationship_ln;
    summary.trio_me = selected.trio_me;
    summary.trio_me_pass = selected.trio_me_pass ? 1 : 0;
    summary.near_tie = near_tie ? 1 : 0;
    summary.rejection_reason = rejection_reason;
    summary.graph_revision_before = graph_before;
    summary.graph_revision_after = graph_before;
    summary.dam_candidate_count = dam_count;
    summary.sire_candidate_count = sire_count;
    summary.configurations_total = static_cast<int>(configurations.size());
    summary.single_configurations_total = single_total;
    summary.pair_configurations_total = pair_total;
    summary.single_hard_feasible = single_hard;
    summary.pair_hard_feasible = pair_hard;
    summary.single_trio_me_pass = single_me;
    summary.pair_trio_me_pass = pair_me;
    summary.single_relationship_pass = single_relationship;
    summary.pair_relationship_pass = pair_relationship;
    summary.single_eligible = single_eligible;
    summary.pair_eligible = pair_eligible;
    summary.configurations_eligible = single_eligible + pair_eligible + 1;
    summary.current_hard_valid = current.hard_valid ? 1 : 0;
    summary.current_assignment_confident =
      current.assignment_confident ? 1 : 0;
    summary.current_retention_eligible =
      current.eligible_for_retention ? 1 : 0;
    summary.recheck_evaluation = input.recheck_only ? 1 : 0;

    if(input.private_frozen_decision_audit){
      // Counterfactual acceptance only: do not change the graph, epochs,
      // messages, candidate overlays or accepted-move callbacks. Every visit
      // sees the SAME supplied final state, including after a would-accept.
      summary.accepted = rejection_reason==0 &&
        (selected.dam!=starting_dam || selected.sire!=starting_sire);
    } else if (rejection_reason == 0 &&
        (selected.dam != starting_dam || selected.sire != starting_sire)) {
      AcceptedMove move;
      move.child = child;
      move.old_dam = starting_dam;
      move.old_sire = starting_sire;
      move.new_dam = selected.dam;
      move.new_sire = selected.sire;
      move.current_objective_ln = input.joint_identity_objective ? decision_objective(current) : current.objective_ln;
      move.selected_objective_ln = input.joint_identity_objective ? decision_objective(selected) : selected.objective_ln;
      move.objective_delta_ln = configuration_contrast(selected,current,input.joint_identity_objective);
      move.dam_support_ln = selected.dam_support_ln;
      move.sire_support_ln = selected.sire_support_ln;
      move.pair_support_ln = selected.pair_support_ln;
      move.identity_margin_ln = selected.identity_margin_ln;
      move.trio_me = selected.trio_me;
      move.graph_revision_before = graph_before;

      // Do not retain a structural plan through a graph patch.  Message
      // revision remains unchanged here by the source-faithful pass contract.

      const DirtyRegion dirty = state.apply_move(
        child, selected.dam, selected.sire);
      if(pap_internal::native_demand::live)pap_internal::native_demand::live->accepted_move(child,starting_dam,starting_sire,selected.dam,selected.sire,state.pass.dam,state.pass.sire,state.graph_revision);
      const PedigreeGraphSnapshot graph = state.graph_snapshot();
      {
        worker.synchronize_graph(graph, dirty);
      }

      move.graph_revision_after = state.graph_revision;
      summary.graph_revision_after = state.graph_revision;
      summary.accepted = 1;
      moves.push_back(move);

    }
    summary.focal_wall_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - focal_started).count();

    summary.message_scoring_model = input.message_scoring ?
      (input.stored_message_scoring ? 2 : 1) : 0;
    for (const auto& configuration : configurations)
      if (configuration.raw_message_approximate) summary.raw_message_approximate = 1;
    summaries.push_back(summary);
    // A focal exclusion is part of plan identity: rebuild even if the graph
    // revision stayed fixed.

    }
  }
  report_visits(input, focal_window_size, focal_window_size, moves.size());
  input.message_scoring = requested_message_model;

  {

    merge_visit_counters(compiled.telemetry, worker.pass.telemetry);
  }

  const double elapsed = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - started).count();

  const auto result_materialization_started = std::chrono::steady_clock::now();

  IntegerMatrix final_parents(input.n, 2);
  IntegerVector final_role(input.n);
  for (int node = 1; node <= input.n; ++node) {
    final_parents(node - 1, 0) = state.pass.dam[node];
    final_parents(node - 1, 1) = state.pass.sire[node];
    // R records an inferred role only for unknown-sex individuals; known-sex
    // role consistency is validated but does not populate this vector.
    final_role[node - 1] = is_na(input.sex[node - 1]) ?
      (state.role_counts[node][0] > 0 ? 1 :
       (state.role_counts[node][1] > 0 ? 2 : 0)) : 0;
  }
  DataFrame accepted_moves = accepted_move_rows(moves);
  DataFrame focal_summary = focal_summary_rows(summaries);
  List output = List::create(
    _["schema_version"] = "simp-parentage-pass-v1",
    _["final_parents"] = final_parents,
    _["final_inferred_role"] = final_role,
    _["accepted_moves"] = accepted_moves,
    _["focal_summary"] = focal_summary,
    _["final_graph_revision"] = state.graph_revision,
    _["message_revision"] = input.message_revision,
    _["locus_parallel_provenance"] = worker_team.provenance());
  if(!input.diagnostic_visit_order.empty()) {
    output["diagnostic_visit_rows"] = visit_diagnostic_rows(diagnostic_visit_rows);
    output["diagnostic_scored_panels"] = panel_diagnostic_rows(diagnostic_panels);
  }
  output["typed_age_changed_revisions"]=wrap(state.age_changed_revision);
  output["typed_age_propagation"]=NumericVector::create(
    _["proposals"]=state.age_proposals,_["rejected_proposals"]=state.age_rejections,
    _["proposal_seconds"]=state.age_proposal_seconds,_["moves"]=state.age_moves,
    _["move_seconds"]=state.age_move_seconds,_["move_component_nodes"]=state.age_move_nodes,
    _["move_relaxations"]=state.age_move_relaxations,_["visit_seconds"]=state.age_visit_seconds);

  std::vector<int> dirty_union_ids;
  for (std::size_t node = 1;
       node < state.publication_dirty_union.size(); ++node)
    if (state.publication_dirty_union[node])
      dirty_union_ids.push_back(static_cast<int>(node));
  output["dirty_union_node_ids"] = wrap(dirty_union_ids);

  output["search_loop_engine"] = input.search_loop_engine;

  if (input.configuration_diagnostics) {
    DataFrame configuration_frame = configuration_rows(configuration_diagnostics, input);
    output["configuration_diagnostics"] = configuration_frame;

    DataFrame provenance_frame = relationship_provenance_rows(configuration_diagnostics);
    output["relationship_provenance"] = provenance_frame;

    output["diagnostic_graph_states"] = graph_diagnostic_rows(graph_state_diagnostics, input.n);
  }
  if (input.configuration_locus_diagnostics) {
    output["configuration_locus_diagnostics"] = locus_diagnostic_rows(configuration_locus_diagnostics);
  }

  return output;
}

}
