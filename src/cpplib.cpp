#include<pybind11/pybind11.h>
#include<pybind11/stl.h>
#include<pybind11/chrono.h>
#include<pybind11/complex.h>
#include<pybind11/functional.h>

#include "Bio.hpp"
#include "include/pybind11_json.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cpplib, m) {
    m.doc() = "c++ extension to python code"; // optional module docstring
    py::class_<Bio>(m, "Bio")
    .def(py::init<>())
    .def("toCodon", &Bio::toCodons)
    .def("reverseComplement", &Bio::reverse_comp)
    .def("ORFfinder", &Bio::ORF_finder, py::arg("gene"), py::arg("minNuc")=70)
    .def("codonCounter", &Bio::codonCounter)
    .def("aminoAcidsCounter", &Bio::aminoAcidsCounter)
    .def("translate", &Bio::translate)
    .def("nucCounter", &Bio::nucleotideFrequency)
    .def_readwrite("codons", &Bio::codons)
    .def_readwrite("gene", &Bio::gene)
    ;
}
