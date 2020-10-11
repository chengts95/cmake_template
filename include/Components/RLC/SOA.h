#pragma once
#include <libflatarray/flat_array.hpp>
#include "RLC.h"



LIBFLATARRAY_REGISTER_SOA(
    Capacitor,
    ((Real)(C))
    ((int2)(port))
    ((Real)(Ieq))
    ((Real)(Geq))
  
    )