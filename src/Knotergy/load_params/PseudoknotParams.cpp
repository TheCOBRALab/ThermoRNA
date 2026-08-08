#include "PseudoknotParams.hpp"
#include "Knotergy/utils/FileUtils.hpp"
#include "Knotergy/utils/common.hpp"
#include "Knotergy/utils/colors.hpp"

namespace thermorna::knotergy {

pk_param PseudoknotParams::load_pk_param(const std::string& paramFile, RoundMethod round_method) {
    ParamSourceInfo info;
    info.label = "Pseudoknot";
    info.requested_path = paramFile;

    if (paramFile.empty()) {
        std::cout << WARNING << " Warning: Default pseudoknot parameter file not found. Using "
        "hard-coded defaults.\n";

        pk_param pkp = pk_param();
        info.status = ParamStatus::Defaulted;
        pkp.set_source_info(info);  // Set the source info for reporting purposes
        return pkp;
    }

    if (!FileUtils::file_exists(paramFile)) {
        THROW_ERROR("Pseudoknot parameters JSON file \"" + paramFile + "\" not found.");
    }

    if (!FileUtils::is_file(paramFile)) {
        THROW_ERROR("Pseudoknot parameters path \"" + paramFile + "\" is not a file.");
    }

    pk_param pkp = parse_pk_json(paramFile, round_method);

    info.status = ParamStatus::LoadedUserFile;
    info.resolved_path = paramFile;
    info.resolved_name = pkp.name;

    pkp.set_source_info(info);

    return pkp;
}

pk_param PseudoknotParams::parse_pk_json(const std::string& jsonFile, RoundMethod round_method) {
    std::ifstream f(jsonFile);
    if (!f.is_open()) {
        THROW_ERROR("Error: Unable to open pseudoknot parameter file: " + jsonFile);
    }

    json data = json::parse(f);

    auto it = data.find("pseudoknot_parameters");
    if (it == data.end()) {
        THROW_ERROR("Missing required key 'pseudoknot_parameters' in file " + jsonFile);
    }
    const auto& pk = *it;

    return pk_param(pk.value("name", std::string{"No Name Provided"}), pk.at("pk_in_ext").get<int>(),
                    pk.at("pk_in_mloop").get<int>(), pk.at("pk_in_pk").get<int>(),
                    pk.at("band_penalty").get<int>(), pk.at("unpaired_in_pk").get<int>(),
                    pk.at("cr_in_pk").get<int>(), pk.at("pk_stack_x").get<double>(),
                    pk.at("pk_internal_x").get<double>(), pk.at("pk_mloop_init").get<int>(),
                    pk.at("pk_mloop_bp").get<int>(), pk.at("pk_mloop_unpaired").get<int>(),
                    round_method);
    }
}