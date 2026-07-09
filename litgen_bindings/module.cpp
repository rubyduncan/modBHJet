#include <nanobind/nanobind.h>

namespace nb = nanobind;

void py_init_module_modbhjet(nb::module_ &m);

NB_MODULE(bhjet, m)
{
    py_init_module_modbhjet(m);
}
