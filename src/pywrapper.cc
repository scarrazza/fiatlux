// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

// namespaces
#include "fiatlux/fiatlux.h"

// namespaces
using namespace fiatlux;
namespace py = pybind11;

// macros
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

PYBIND11_MODULE(_core, m) {
    py::class_<FiatLux>(m, "FiatLux")
        .def(py::init<const std::string &>())
        .def("PlugAlphaQED", &FiatLux::PlugAlphaQED)
        .def("PlugStructureFunctions", &FiatLux::PlugStructureFunctions)
        .def("InsertInelasticSplitQ", &FiatLux::InsertInelasticSplitQ)
        .def("EvaluatePhoton", &FiatLux::EvaluatePhoton);

    py::class_<luxqed>(m, "luxqed")
        .def(py::init())
        .def_readonly("elastic", &luxqed::elastic)
        .def_readonly("inelastic_pf", &luxqed::inelastic_pf)
        .def_readonly("msbar_pf", &luxqed::msbar_pf)
        .def_readonly("total", &luxqed::total);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
