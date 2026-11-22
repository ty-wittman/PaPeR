#include "search_inputs.h"
namespace pap_internal::search {

auto SearchInputs::elide_raw() const -> bool  {
    return skip_unused_raw && identity_prescreen && !identity_union &&
      pap_internal::native_demand::live && config->likelihood_ancestral &&
      joint_identity_objective && joint_support_gates && null_joint_identity &&
      relationship_mode == 1 && lazy_independent_gates &&
      config->evaluation.relationship_common_scaled_fast && gate_window_ln == 0.0 &&
      (!configuration_diagnostics || (private_frozen_decision_audit && private_decision_only_diagnostics)) &&
      !configuration_locus_diagnostics &&
      !raw_joint_census && !raw_support_observer &&
      !benchmark_baseline_mode &&
      !message_scoring &&
      !private_test_raw_comparison_exception && !private_test_relationship_comparison_exception &&
      !private_test_configuration_post_score_exception;
  }

auto SearchInputs::int_field(const List& inputs, const char* name) -> IntegerVector  {
    if (!inputs.containsElementNamed(name))
      engine_stop("whole-pass input is missing '%s'", name);
    return as<IntegerVector>(inputs[name]);
  }

auto SearchInputs::numeric_field(const List& inputs, const char* name) -> NumericVector  {
    if (!inputs.containsElementNamed(name))
      engine_stop("whole-pass input is missing '%s'", name);
    return as<NumericVector>(inputs[name]);
  }

SearchInputs::SearchInputs(const List& inputs, int n_) : config(pap_internal::native_demand::live?pap_internal::native_demand::live->execution_config():nullptr), n(n_) {
    if(!config)config=std::make_shared<const pap_internal::ExecutionConfig>(inputs);
    full_contextual_relationship_index.fill(-1);
    const IntegerVector focal = int_field(inputs, "focal_order");
    focal_order.assign(focal.begin(), focal.end());
    if(inputs.containsElementNamed("diagnostic_visit_order"))
      diagnostic_visit_order=as<std::vector<int>>(inputs["diagnostic_visit_order"]);
    if(inputs.containsElementNamed("diagnostic_visit_recheck_order"))
      diagnostic_visit_recheck_order=as<std::vector<int>>(inputs["diagnostic_visit_recheck_order"]);
    const IntegerVector d_off = int_field(inputs, "dam_offsets");
    const IntegerVector s_off = int_field(inputs, "sire_offsets");
    const IntegerVector d_ids = int_field(inputs, "dam_ids");
    const IntegerVector s_ids = int_field(inputs, "sire_ids");
    const NumericVector d_quick = numeric_field(inputs, "dam_quick");
    const NumericVector s_quick = numeric_field(inputs, "sire_quick");
    const IntegerVector rank = int_field(inputs, "id_order_rank");
    const IntegerVector sex_in = int_field(inputs, "sex");
    const IntegerVector lo = int_field(inputs, "birth_lo");
    const IntegerVector hi = int_field(inputs, "birth_hi");
    const IntegerVector last = int_field(inputs, "year_last");
    const NumericVector min_age_in = numeric_field(inputs, "min_age");
    const NumericVector max_age_in = numeric_field(inputs, "max_age");
    dam_offsets.assign(d_off.begin(), d_off.end());
    sire_offsets.assign(s_off.begin(), s_off.end());
    dam_ids.assign(d_ids.begin(), d_ids.end());
    sire_ids.assign(s_ids.begin(), s_ids.end());
    dam_quick.assign(d_quick.begin(), d_quick.end());
    sire_quick.assign(s_quick.begin(), s_quick.end());
    id_order_rank.assign(rank.begin(), rank.end());
    sex.assign(sex_in.begin(), sex_in.end());
    birth_lo.assign(lo.begin(), lo.end());
    birth_hi.assign(hi.begin(), hi.end());
    year_last.assign(last.begin(), last.end());
    min_age.assign(min_age_in.begin(), min_age_in.end());
    max_age.assign(max_age_in.begin(), max_age_in.end());
    Tassign_ln = as<double>(inputs["Tassign_ln"]);
    identity_margin_ln = as<double>(inputs["identity_margin_ln"]);
    tol = as<double>(inputs["tol"]);
    if (inputs.containsElementNamed("tie_tolerance"))
      tie_tolerance = as<double>(inputs["tie_tolerance"]);
    max_mendelian_errors = as<int>(inputs["max_mendelian_errors"]);
    const bool ancestral_prescreen_default=inputs.containsElementNamed("pedigree_likelihood")&&as<bool>(inputs["pedigree_likelihood"])&&
      (!inputs.containsElementNamed("likelihood_evidence")||as<std::string>(inputs["likelihood_evidence"])=="ancestral");
    const bool prescreen_shadow=inputs.containsElementNamed("pair_rescue_shadow")&&as<bool>(inputs["pair_rescue_shadow"]);
    const bool prescreen_calibration=inputs.containsElementNamed("pair_rescue_calibration")&&as<bool>(inputs["pair_rescue_calibration"]);
    if(ancestral_prescreen_default){
      pair_rescue_enabled=pap_internal::native_demand::live&&!prescreen_shadow&&!prescreen_calibration;
      pair_rescue_support_only=true;pair_rescue_margin_ln=std::log(100.0);pair_rescue_top_k=10;
    }
    if(inputs.containsElementNamed("pair_rescue"))pair_rescue_enabled=as<bool>(inputs["pair_rescue"]);
    if(inputs.containsElementNamed("complete_identity_competitors"))complete_identity_competitors=as<bool>(inputs["complete_identity_competitors"]);
    if(inputs.containsElementNamed("completion_rank_only"))completion_rank_only=as<bool>(inputs["completion_rank_only"]);
    if(inputs.containsElementNamed("likelihood_raw_joint_census"))raw_joint_census=as<bool>(inputs["likelihood_raw_joint_census"]);
    if(inputs.containsElementNamed("pair_rescue_shadow"))pair_rescue_shadow=as<bool>(inputs["pair_rescue_shadow"]);
    if(inputs.containsElementNamed("pair_rescue_diagnostics"))pair_rescue_diagnostics=as<bool>(inputs["pair_rescue_diagnostics"]);
    if(inputs.containsElementNamed("pair_rescue_support_only"))pair_rescue_support_only=as<bool>(inputs["pair_rescue_support_only"]);
    if(inputs.containsElementNamed("pair_rescue_calibration"))pair_rescue_calibration=as<bool>(inputs["pair_rescue_calibration"]);
    if(pair_rescue_calibration&&(!pair_rescue_shadow||!pair_rescue_support_only))
      engine_stop("pair-rescue calibration requires support-only shadow mode");
    if(inputs.containsElementNamed("pair_rescue_margin_ln"))pair_rescue_margin_ln=as<double>(inputs["pair_rescue_margin_ln"]);
    if(inputs.containsElementNamed("pair_rescue_top_k")){
      const double k=as<double>(inputs["pair_rescue_top_k"]);
      if(!std::isfinite(k)||k<1||k>100000||std::floor(k)!=k)engine_stop("stage1_pair_rescue_top_k must be a positive integer");
      pair_rescue_top_k=int(k);
    }
    if(std::isnan(pair_rescue_margin_ln)||pair_rescue_margin_ln<0)engine_stop("stage1_pair_rescue_margin_ln must be nonnegative");
    if(pair_rescue_enabled||pair_rescue_shadow){
      if(!inputs.containsElementNamed("pedigree_likelihood")||!as<bool>(inputs["pedigree_likelihood"])||!pap_internal::native_demand::live)
        engine_stop("candidate pair rescue requires the M_ped2 live scorer");
      if(pair_rescue_enabled&&pair_rescue_shadow)engine_stop("pair rescue and its shadow are separate modes");
      if(!inputs.containsElementNamed("Tfilter_ln")||(std::isnan(as<double>(inputs["Tfilter_ln"]))||as<double>(inputs["Tfilter_ln"])==R_PosInf))
        engine_stop("pair rescue requires Tfilter_ln below positive infinity for its trio screen");
      pair_rescue_trio_ln=2*as<double>(inputs["Tfilter_ln"]);
    }
    if (inputs.containsElementNamed("trio_lr_margin_ln")) {
      const double value = as<double>(inputs["trio_lr_margin_ln"]);
      trio_lr_margin_ln = ISNAN(value) ? 0.0 : value;
      if (trio_lr_margin_ln < 0.0) engine_stop("stage1_trio_lr_margin_ln must be nonnegative");
    }
    if (inputs.containsElementNamed("gate_window_full"))
      gate_window_full = as<bool>(inputs["gate_window_full"]);
    if (inputs.containsElementNamed("gate_window_ln")) {
      const double value = as<double>(inputs["gate_window_ln"]);
      gate_window_ln = ISNAN(value) ? 0.0 : value;
      if (gate_window_ln < 0.0) engine_stop("stage1_gate_window_ln must be nonnegative");
      gate_window_base = gate_window_ln;
    }
    if (inputs.containsElementNamed("duo_lr_margin_ln")) {
      const double value = as<double>(inputs["duo_lr_margin_ln"]);
      duo_lr_margin_ln = ISNAN(value) ? 0.0 : value;
      if (duo_lr_margin_ln < 0.0) engine_stop("stage1_duo_lr_margin_ln must be nonnegative");
    }
    if (inputs.containsElementNamed("duo_lr_margin_recheck_only"))
      duo_lr_margin_recheck_only = as<bool>(inputs["duo_lr_margin_recheck_only"]);
    if (inputs.containsElementNamed("trio_lr_margin_recheck_only"))
      trio_lr_margin_recheck_only = as<bool>(inputs["trio_lr_margin_recheck_only"]);
    if (inputs.containsElementNamed("trio_lr_filter_ln")) {
      const double value = as<double>(inputs["trio_lr_filter_ln"]);
      trio_lr_filter_ln = ISNAN(value) ? R_NegInf : value;
    }
    start_graph_revision = as<int>(inputs["start_graph_revision"]);
    message_revision = as<int>(inputs["message_revision"]);
    if (inputs.containsElementNamed("configuration_diagnostics"))
      configuration_diagnostics = as<bool>(
        inputs["configuration_diagnostics"]);
    if (inputs.containsElementNamed("configuration_locus_diagnostics"))
      configuration_locus_diagnostics = as<bool>(
        inputs["configuration_locus_diagnostics"]);
    if (inputs.containsElementNamed("diagnostic_start_position"))
      diagnostic_start_position = as<int>(inputs["diagnostic_start_position"]);
    if (inputs.containsElementNamed("diagnostic_end_position"))
      diagnostic_end_position = as<int>(inputs["diagnostic_end_position"]);
    if (inputs.containsElementNamed("initial_dirty_union_node_ids"))
      initial_dirty_union_node_ids = as<std::vector<int>>(
        inputs["initial_dirty_union_node_ids"]);
    if (inputs.containsElementNamed("recheck_only"))
      recheck_only = as<bool>(inputs["recheck_only"]);
    if(inputs.containsElementNamed("private_decision_only_diagnostics"))
      private_decision_only_diagnostics=as<bool>(inputs["private_decision_only_diagnostics"]);
    if(inputs.containsElementNamed("private_frozen_decision_audit"))
      private_frozen_decision_audit=as<bool>(inputs["private_frozen_decision_audit"]);
    // Summary-only frozen audits preserve lazy gates; requesting full
    // configuration diagnostics deliberately resolves extra comparisons.
    if (inputs.containsElementNamed("recheck_focals")) {
      const IntegerVector focals = as<IntegerVector>(inputs["recheck_focals"]);
      recheck_focals.assign(focals.begin(), focals.end());
      for (int focal : recheck_focals)
        if (focal < 1 || focal > n)
          engine_stop("complete recheck focal is outside the native ID range");
    }
    if (inputs.containsElementNamed("configuration_locus_focal_position"))
      configuration_locus_focal_position = as<int>(
        inputs["configuration_locus_focal_position"]);

    if (inputs.containsElementNamed("private_test_relationship_comparison_exception"))
      private_test_relationship_comparison_exception = as<bool>(
        inputs["private_test_relationship_comparison_exception"]);
    if (inputs.containsElementNamed("private_complete_message_trial"))
      private_complete_message_trial = as<bool>(inputs["private_complete_message_trial"]);

    if (inputs.containsElementNamed("private_test_raw_comparison_exception"))
      private_test_raw_comparison_exception = as<bool>(
        inputs["private_test_raw_comparison_exception"]);
    if (inputs.containsElementNamed("private_test_configuration_post_score_exception"))
      private_test_configuration_post_score_exception = as<bool>(
        inputs["private_test_configuration_post_score_exception"]);
    if (inputs.containsElementNamed("lazy_independent_gates"))
      lazy_independent_gates = as<bool>(inputs["lazy_independent_gates"]);
    if (inputs.containsElementNamed("incumbent_full_scoring"))
      incumbent_full_scoring = as<bool>(inputs["incumbent_full_scoring"]);
    joint_identity_objective=pap_internal::native_demand::live!=nullptr;
    if(inputs.containsElementNamed("identity_objective")) {
      const std::string mode=as<std::string>(inputs["identity_objective"]);
      if(mode!="joint"&&mode!="legacy_side_u")engine_stop("stage1_identity_objective must be joint or legacy_side_u");
      joint_identity_objective=mode=="joint";
    }
    if(inputs.containsElementNamed("joint_support_gates"))joint_support_gates=as<bool>(inputs["joint_support_gates"]);
    if(inputs.containsElementNamed("prescreen_identity"))identity_prescreen=as<bool>(inputs["prescreen_identity"]);
    if(inputs.containsElementNamed("likelihood_skip_unused_raw"))skip_unused_raw=as<bool>(inputs["likelihood_skip_unused_raw"]);
    if(inputs.containsElementNamed("likelihood_identity_union"))identity_union=as<bool>(inputs["likelihood_identity_union"]);
    if(inputs.containsElementNamed("pedigree_joint_support_diagnostics"))raw_support_observer=as<bool>(inputs["pedigree_joint_support_diagnostics"]);
    if(inputs.containsElementNamed("support_available_baselines"))support_available_baselines=as<bool>(inputs["support_available_baselines"]);
    if(inputs.containsElementNamed("null_joint_identity"))null_joint_identity=as<bool>(inputs["null_joint_identity"]);
    // I(empty)=0 is a property of configuration-local likelihood ratios.
    // The retained union/full references use absolute joint log likelihoods;
    // their null must stay on that same scale as every competing row.
    null_joint_identity = null_joint_identity && config->likelihood_ancestral && !identity_union;
    if (inputs.containsElementNamed("raw_common_background_legacy_missing_mates"))
      raw_common_background_legacy_missing_mates = as<bool>(
        inputs["raw_common_background_legacy_missing_mates"]);
    if (inputs.containsElementNamed("raw_common_observation_offsets")) {
      raw_common_observation_offsets = as<std::vector<double>>(
        inputs["raw_common_observation_offsets"]);
      if (raw_common_observation_offsets.size() != static_cast<std::size_t>(n) ||
          !std::all_of(raw_common_observation_offsets.begin(),
            raw_common_observation_offsets.end(), [](double x) { return std::isfinite(x); }))
        engine_stop("raw common-background observation offsets must cover all typed IDs");
    }
    if (inputs.containsElementNamed("message_scoring"))
      message_scoring = as<bool>(inputs["message_scoring"]);
    if (inputs.containsElementNamed("stored_message_scoring"))
      stored_message_scoring = as<bool>(inputs["stored_message_scoring"]);
    message_scoring = message_scoring || stored_message_scoring;
    if (inputs.containsElementNamed("raw_scaled_fast"))
      raw_scaled_fast = as<bool>(inputs["raw_scaled_fast"]);
    if (inputs.containsElementNamed("adaptive_selective_publisher"))
      adaptive_selective_publisher = as<bool>(
        inputs["adaptive_selective_publisher"]);
    if (inputs.containsElementNamed("adaptive_selective_threshold"))
      adaptive_selective_threshold = as<double>(
        inputs["adaptive_selective_threshold"]);
    if (inputs.containsElementNamed("progress_parentage_round"))
      progress_parentage_round = as<int>(inputs["progress_parentage_round"]);
    if (inputs.containsElementNamed("progress_locus_interval"))
      progress_locus_interval = as<int>(inputs["progress_locus_interval"]);
    if (inputs.containsElementNamed("progress_heartbeat_seconds"))
      progress_heartbeat_seconds = as<double>(
        inputs["progress_heartbeat_seconds"]);
    if (inputs.containsElementNamed("progress_locus_heartbeats"))
      progress_locus_heartbeats = as<bool>(
        inputs["progress_locus_heartbeats"]);
    if (inputs.containsElementNamed("progress_callback"))
      progress_callback = inputs["progress_callback"];
    if (inputs.containsElementNamed("visit_progress_callback"))
      visit_progress_callback = inputs["visit_progress_callback"];
    if (inputs.containsElementNamed("benchmark_baseline_mode"))
      benchmark_baseline_mode = as<bool>(inputs["benchmark_baseline_mode"]);
    if (inputs.containsElementNamed("benchmark_opt5_diag_mode"))
      benchmark_opt5_diag_mode = as<bool>(inputs["benchmark_opt5_diag_mode"]);
    if (inputs.containsElementNamed("phase_timing"))
      phase_timing = as<bool>(inputs["phase_timing"]);

    if (inputs.containsElementNamed("search_loop_engine"))
      search_loop_engine = as<std::string>(inputs["search_loop_engine"]);
    if (inputs.containsElementNamed("relationship_mode"))
      relationship_mode = as<int>(inputs["relationship_mode"]);
    if (inputs.containsElementNamed("full_relationship_registry")) {
      const List registry = as<List>(inputs["full_relationship_registry"]);
      const CharacterVector relationship_id = as<CharacterVector>(registry["relationship_id"]);
      const NumericMatrix delta = as<NumericMatrix>(registry["delta"]);
      const IntegerVector orientation_mask = as<IntegerVector>(
        registry["orientation_mask"]);
      const LogicalVector parent = as<LogicalVector>(registry["parent_numerator"]);
      const LogicalVector alternative = as<LogicalVector>(
        registry["denominator_alternative"]);
      if (delta.ncol() != 9 || delta.nrow() != relationship_id.size() ||
          delta.nrow() != orientation_mask.size() ||
          delta.nrow() != parent.size() ||
          delta.nrow() != alternative.size())
        engine_stop("native full relationship registry payload is malformed");
      full_relationship_delta.resize(static_cast<std::size_t>(delta.nrow()));
      full_relationship_orientation_mask.resize(static_cast<std::size_t>(delta.nrow()));
      full_relationship_parent_numerator.resize(static_cast<std::size_t>(delta.nrow()));
      full_relationship_denominator_alternative.resize(static_cast<std::size_t>(delta.nrow()));
      for (int row = 0; row < delta.nrow(); ++row) {
        const std::string id = as<std::string>(relationship_id[row]);
        if (id.empty() || full_relationship_index.count(id))
          engine_stop("native full relationship registry IDs are invalid");
        full_relationship_index[id] = row;
        for (int column = 0; column < 9; ++column)
          full_relationship_delta[static_cast<std::size_t>(row)][column] = delta(row, column);
        if (orientation_mask[row] < 0 || orientation_mask[row] > 3)
          engine_stop("native full relationship orientation mask is invalid");
        full_relationship_orientation_mask[static_cast<std::size_t>(row)] =
          static_cast<unsigned char>(orientation_mask[row]);
        full_relationship_parent_numerator[static_cast<std::size_t>(row)] = parent[row];
        full_relationship_denominator_alternative[static_cast<std::size_t>(row)] = alternative[row];
      }
      static constexpr std::array<const char*, 9> contextual_id{{
        "HS_HA", "HS_F1C", "HS_FA", "HS_GP", "HS_GGG", "FA_HA",
        "FS_HA", "FS_H1C", "FS_PARENTS_FS"}};
      for (std::size_t i = 0; i < contextual_id.size(); ++i) {
        const auto found = full_relationship_index.find(contextual_id[i]);
        if (found != full_relationship_index.end())
          full_contextual_relationship_index[i] = found->second;
      }
    }
    if (inputs.containsElementNamed("bounded_kernel_profile"))
      bounded_kernel_profile = as<bool>(inputs["bounded_kernel_profile"]);

    if (inputs.containsElementNamed("focal_start_position"))
      focal_start_position = as<int>(inputs["focal_start_position"]);
    if (inputs.containsElementNamed("focal_end_position"))
      focal_end_position = as<int>(inputs["focal_end_position"]);
    if (inputs.containsElementNamed("replay_moves"))
      replay_moves = as<IntegerMatrix>(inputs["replay_moves"]);
    if (inputs.containsElementNamed("shortlist_mode")) {
      const std::string mode = as<std::string>(inputs["shortlist_mode"]);
      fixed50 = mode == "sequoia50";
    }
  }

auto SearchInputs::is_recheck_focal(int focal) const -> bool  {
    return !recheck_only || std::find(recheck_focals.begin(),
                                      recheck_focals.end(), focal) !=
      recheck_focals.end();
  }

auto SearchInputs::validate() const -> void  {
    if (n < 1) engine_stop("whole-pass input requires at least one individual");

    if (relationship_mode != 0 && relationship_mode != 1)
      engine_stop("whole-pass relationship_mode must be 0 (simp) or 1 (full)");
    if (relationship_mode == 1 && full_relationship_delta.empty())
      engine_stop("whole-pass full relationship mode requires a registry payload");
    if (!full_relationship_delta.empty()) for (const char* id : {
        "U", "PO", "FS", "HS", "GP", "FA", "HA", "HS_PO", "PO_HA",
        "HS_HA", "HS_F1C", "HS_FA", "HS_GP", "HS_GGG", "FA_HA", "FS_HA",
        "FS_H1C", "FS_PARENTS_FS"}) if (!full_relationship_index.count(id))
      engine_stop("native full relationship registry is missing a required identity");
    if (relationship_mode == 1) for (int index : full_contextual_relationship_index)
      if (index < 0 || index >= static_cast<int>(full_relationship_delta.size()) ||
          full_relationship_orientation_mask[static_cast<std::size_t>(index)] == 0)
        engine_stop("native full contextual relationship registry is malformed");
    if (static_cast<int>(focal_order.size()) != n)
      engine_stop("whole-pass focal order must contain every individual exactly once");
    if (static_cast<int>(dam_offsets.size()) != n + 1 ||
        static_cast<int>(sire_offsets.size()) != n + 1)
      engine_stop("whole-pass candidate offsets must have length n + 1");
    if (dam_offsets.front() != 0 || sire_offsets.front() != 0 ||
        dam_offsets.back() != static_cast<int>(dam_ids.size()) ||
        sire_offsets.back() != static_cast<int>(sire_ids.size()))
      engine_stop("whole-pass candidate offsets are inconsistent with flat arrays");
    auto validate_offsets = [&](const std::vector<int>& offsets,
                                const std::vector<int>& ids,
                                const std::vector<double>& quick,
                                const char* label) {
      for (int i = 0; i < n; ++i) {
        if (offsets[i] < 0 || offsets[i] > offsets[i + 1] ||
            offsets[i + 1] > static_cast<int>(ids.size()))
          engine_stop("whole-pass %s offsets are not monotone", label);
      }
      if (quick.size() != ids.size())
        engine_stop("whole-pass %s candidate metadata is not aligned", label);
      for (int id : ids) if (id < 1 || id > n)
        engine_stop("whole-pass %s contains an out-of-range ID", label);
      for (int focal = 0; focal < n; ++focal) {
        std::vector<int> seen;
        for (int i = offsets[focal]; i < offsets[focal + 1]; ++i) {
          if (std::find(seen.begin(), seen.end(), ids[i]) != seen.end())
            engine_stop("whole-pass %s contains a duplicate candidate", label);
          seen.push_back(ids[i]);
        }
        // The retained universe is fixed at 50, with one explicit exception:
        // a currently assigned parent may be forced back into the compiled
        // pass when it falls just outside that retained top-50 list.
        if (fixed50 && offsets[focal + 1] - offsets[focal] > 51)
          engine_stop("fixed-50 whole-pass candidate set exceeds 50 plus forced current parent for %s",
               label);
      }
    };
    validate_offsets(dam_offsets, dam_ids, dam_quick, "dam");
    validate_offsets(sire_offsets, sire_ids, sire_quick, "sire");
    auto validate_permutation = [&](const std::vector<int>& values,
                                    const char* label) {
      if (static_cast<int>(values.size()) != n)
        engine_stop("whole-pass %s must have length n", label);
      std::vector<int> seen(static_cast<std::size_t>(n + 1), 0);
      for (int value : values) {
        if (value < 1 || value > n || ++seen[value] != 1)
          engine_stop("whole-pass %s must be a permutation of 1..n", label);
      }
    };
    validate_permutation(focal_order, "focal order");
    if(diagnostic_visit_order.size()>512U)
      engine_stop("private diagnostic visit order exceeds512visits");
    for(int id:diagnostic_visit_order)if(id<1 || id>n)
      engine_stop("private diagnostic visit order contains an out-of-range ID");
    if(!diagnostic_visit_recheck_order.empty() && diagnostic_visit_recheck_order.size()!=diagnostic_visit_order.size())
      engine_stop("private diagnostic recheck order must align with visit order");
    for(int mode:diagnostic_visit_recheck_order)if(mode!=0 && mode!=1)
      engine_stop("private diagnostic recheck order must contain only0or1");
    validate_permutation(id_order_rank, "ID order rank");
    const int effective_end = focal_end_position == 0 ? n : focal_end_position;
    if (focal_start_position < 1 || focal_start_position > n ||
        effective_end < focal_start_position || effective_end > n)
      engine_stop("whole-pass focal window is outside the focal order");
    if (replay_moves.nrow() > 0 && replay_moves.ncol() != 7)
      engine_stop("whole-pass replay moves must have seven integer columns");
    if (sex.size() != static_cast<std::size_t>(n) ||
        birth_lo.size() != static_cast<std::size_t>(n) ||
        birth_hi.size() != static_cast<std::size_t>(n) ||
        year_last.size() != static_cast<std::size_t>(n))
      engine_stop("whole-pass life-history metadata must have length n");
    if (min_age.size() != 2 || max_age.size() != 2 ||
        !std::isfinite(min_age[0]) || !std::isfinite(min_age[1]) ||
        !std::isfinite(max_age[0]) || !std::isfinite(max_age[1]) ||
        min_age[0] < 0 || min_age[1] < 0 || max_age[0] < min_age[0] ||
        max_age[1] < min_age[1])
      engine_stop("whole-pass age bounds are invalid");
    const auto valid_log_threshold = [](double value) {
      return R_finite(value) || is_negative_infinity_ln(value);
    };
    if (!valid_log_threshold(Tassign_ln) ||
        !valid_log_threshold(identity_margin_ln) ||
        !std::isfinite(tol) || tol < 0 || !std::isfinite(tie_tolerance) ||
        tie_tolerance < 0 ||
        progress_locus_interval < 1 ||
        !std::isfinite(progress_heartbeat_seconds) ||
        progress_heartbeat_seconds < 0.0 ||
        max_mendelian_errors < 0 || start_graph_revision < 0 ||
        message_revision < 0 ||
        !std::isfinite(adaptive_selective_threshold) ||
        adaptive_selective_threshold < 0.0 ||
        adaptive_selective_threshold > 1.0)
      engine_stop("whole-pass numerical or revision metadata is invalid");
    if (progress_callback != R_NilValue && !Rf_isFunction(progress_callback))
      engine_stop("whole-pass progress callback must be a function or NULL");
    if (visit_progress_callback != R_NilValue &&
        !Rf_isFunction(visit_progress_callback))
      engine_stop("visit progress callback must be a function or NULL");
    const int diagnostic_end = diagnostic_end_position == 0 ? n :
      diagnostic_end_position;
    if (diagnostic_start_position < 1 || diagnostic_end < diagnostic_start_position ||
        diagnostic_end > n)
      engine_stop("invalid diagnostic focal window");
    for (int node : initial_dirty_union_node_ids)
      if (node < 1 || node > n || !adaptive_selective_publisher)
        engine_stop("invalid initial publication dirty union");
    if (configuration_locus_diagnostics &&
        ((configuration_locus_focal_position < 1 &&
          configuration_locus_focal_position != -1) ||
         configuration_locus_focal_position > n))
      engine_stop("configuration locus diagnostics require a focal position in 1..n or -1 for the diagnostic window");

    const bool bounded_messages = search_loop_engine == "bounded_messages";
    if (search_loop_engine != "visited_cavity" && !bounded_messages)
      engine_stop("simp_loop_engine must be 'visited_cavity' or 'bounded_messages'");

  }

}
