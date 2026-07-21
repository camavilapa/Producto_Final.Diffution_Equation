#pragma once
#include <iostream>

template <class deriv_t, class system_t, class printer_t>
void integrate_euler(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer)
{
    // vector to store derivs
    system_t dsdt(s.size());
    system_t k1(s.size());

    // time loop
    
    for(double t = tinit; t <= tend; t = t + dt) { // NOTE: Last time step not necessarily t

        writer(s, t);
        // compute derivs
        fderiv(s, dsdt, t);
        k1 = dsdt;

        s=s+dt*k1;

      }
}