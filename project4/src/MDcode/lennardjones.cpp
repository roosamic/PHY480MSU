#include "lennardjones.h"
#include "system.h"
#include <cmath>
#include <unitconverter.h>
double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    /*
     This funciton calculates the total force on each atom given by the Lennard Jones potensial
     It also includes a cut-off value rCut where F(r) is more or less zero and an update on the potensial
     energy value
     */

    m_potentialEnergy = 0; // Remember to compute this in the loop
    const int numberOfAtoms = system.atoms().size();
    double epsilon24 = 24*m_epsilon;
    double boxLength = system.systemSize().x();
    double rCut = 2.5*m_sigma;
    double rCut2 = rCut*rCut;
    double sigma6 = m_sigma*m_sigma*m_sigma*m_sigma*m_sigma*m_sigma;
    double sigma12 = sigma6*sigma6;
    double potentialEnergyAtRcut = 4*m_epsilon*(sigma12*pow(rCut,-12) - sigma6*pow(rCut, -6));
    vector<Atom*> &atoms = system.atoms();
    for (int i = 0; i < numberOfAtoms; i++){
        Atom *atom1 = atoms[i];
        double fx = 0;
        double fy = 0;
        double fz = 0;
        double pe = 0;
        double x = atom1->position(0);
        double y = atom1->position(1);
        double z = atom1->position(2);
#pragma omp parallel for reduction(+:fx,fy,fz,pe)
        for (int j = i+1; j < numberOfAtoms; j++){
            Atom *atom2 = atoms[j];

            double dx = x - atom2->position(0);
            double dy = y - atom2->position(1);
            double dz = z - atom2->position(2);

            if ((dx) <= -boxLength/2.0) dx += boxLength;
            if ((dx) > boxLength/2.0) dx -= boxLength;

            if ((dy) <= -boxLength/2.0) dy += boxLength;
            if ((dy) > boxLength/2.0) dy -= boxLength;

            if ((dz) <= -boxLength/2.0) dz += boxLength;
            if ((dz) > boxLength/2.0) dz -= boxLength;

            double dr2 = dx*dx + dy*dy + dz*dz;

            if(dr2 < rCut2) {
                double oneOverDr2 = 1.0 / dr2;
                double oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;
                double oneOverDr12 = oneOverDr6*oneOverDr6;
                double force_scalar = epsilon24*oneOverDr2*(2*sigma12*oneOverDr12 - sigma6*oneOverDr6);

                fx += force_scalar*dx;
                fy += force_scalar*dy;
                fz += force_scalar*dz;

                atom2->force[0] -= force_scalar*dx;
                atom2->force[1] -= force_scalar*dy;
                atom2->force[2] -= force_scalar*dz;

                pe += 4*m_epsilon*(sigma12*oneOverDr12 - sigma6*oneOverDr6) - potentialEnergyAtRcut;
            }
        }
        atom1->force[0] += fx;
        atom1->force[1] += fy;
        atom1->force[2] += fz;
        m_potentialEnergy += pe;
    }
}
