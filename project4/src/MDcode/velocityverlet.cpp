#include "velocityverlet.h"
#include "system.h"
#include "atom.h"

void VelocityVerlet::integrate(System &system, double dt, int N)
{   // Applying the Velocity Verlet algorithm
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*0.5*dt/atom->mass();
        atom->position += atom->velocity*dt;
    }

    system.applyPeriodicBoundaryConditions(N);
    system.calculateForces();

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*0.5*dt/atom->mass();
    }
}
