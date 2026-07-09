from pathlib import Path
import re
import litgen

ROOT = Path(__file__).resolve().parents[1]
print(ROOT)

HEADERS = [
    ROOT / "include" / "JetDynamics.hpp",
    ROOT / "include" / "TargetPhotonField.hpp",
    ROOT / "include" / "TargetBlackBody.hpp",
    ROOT / "include" / "RadiationZone.hpp",
    ROOT / "include" / "BLJet.hpp",
    ROOT / "include" / "IsoJet.hpp",
    ROOT / "include" / "BHJet.hpp",
]

OUT_CPP = ROOT / "litgen_bindings" / "generated_bindings.cpp"
OUT_PYI = ROOT / "_modbhjet.pyi"


def make_options() -> litgen.LitgenOptions:
    options = litgen.LitgenOptions()
    options.bind_library = litgen.BindLibraryType.nanobind
    options.namespaces_root = ["bhjet"]
    options.fn_exclude_by_name__regex = r"reinit_grid_arrays"
    options.member_exclude_by_name__regex = (
        r".*_grid|"
        r"spline_.*|"
        r"jet_dyn|"
        r"nozzle_ener|"
        r"electrons_.*|"
        r"target_component_.*|"
        r"total_target_radiation_.*|"
        r"photon_.*|"
        r"cyclosyn_.*"
    )
    options.member_exclude_by_type__regex = (
        r"gsl_.*|"
        r"kariba::.*|"
        r"jet_dynpars|"
        r"jet_enpars"
    )
    options.class_exclude_by_name__regex = r"jet_dynpars|jet_enpars"
    return options

def patch_generated_cpp() -> None:
        '''for fixing the shared pointers after regenerating the file'''
        text = OUT_CPP.read_text()

        replacements = {
            "#include <nanobind/stl/shared_ptr.h>":
                "#include <nanobind/stl/shared_ptr.h>\n#include <nanobind/stl/pair.h>",
            "nb::class_<bhjet::JetDynamics>":
                "nb::class_<bhjet::JetDynamics>",
            "nb::class_<bhjet::TargetPhotonField>":
                "nb::class_<bhjet::TargetPhotonField>",
            "nb::class_<bhjet::TargetBlackBody>":
                "nb::class_<bhjet::TargetBlackBody, bhjet::TargetPhotonField>",
            "nb::class_<bhjet::BLJet>":
                "nb::class_<bhjet::BLJet, bhjet::JetDynamics>",
            "nb::class_<bhjet::IsoJet>":
                "nb::class_<bhjet::IsoJet, bhjet::JetDynamics>",
            "nb::class_<bhjet::RadiationZone>":
                "nb::class_<bhjet::RadiationZone>",
        }

        for old, new in replacements.items():
            text = text.replace(old, new)

        target_init = '''        .def(nb::init<std::string, double, double, bool, size_t>(),
            nb::arg("name_"), nb::arg("distance_") = defaults::DISTANCE, nb::arg("redshift_") = defaults::REDSHIFT, nb::arg("add_to_total_flux_") = defaults::DEFAULT_ADD_TO_TOTAL_FLUX, nb::arg("verbosity_level_") = defaults::DEFAULT_VERBOSITY_LEVEL,
            " ----------------------------\\n Constructor with defaults\\n ----------------------------")
'''
        text = text.replace(target_init, "")
        text = text.replace("defaults::", "bhjet::defaults::")
        text = re.sub(r'nb::arg\("([A-Za-z][A-Za-z0-9_]*)_"\)', r'nb::arg("\1")', text)

        OUT_CPP.write_text(text)


def main() -> None:
    OUT_CPP.parent.mkdir(parents=True, exist_ok=True)
    OUT_PYI.parent.mkdir(parents=True, exist_ok=True)

    options = make_options()

    litgen.write_generated_code_for_files(
        options=options,
        input_cpp_header_files=[str(h) for h in HEADERS],
        output_cpp_pydef_file=str(OUT_CPP),
        output_stub_pyi_file=str(OUT_PYI),
    )


    patch_generated_cpp()


if __name__ == "__main__":
    main()
