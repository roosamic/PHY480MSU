#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions(int N) {
    // Method for updating the periodiv boundary condition
    // Also updates the initial position of the atoms
    int numberOfAtoms = 4*N*N*N;
    double boxLength = systemSize().x();

    for (int i = 0; i < numberOfAtoms; i++){
        Atom *atom = m_atoms[i];
        for (int dim = 0; dim < 3; dim++){
            if (atom->position(dim) <= 0){
                atom->initial_position(dim) += boxLength;
                atom->position(dim) += boxLength;
            }
            if (atom->position(dim) >= boxLength){
                atom->initial_position(dim) -= boxLength;
                atom->position(dim) -= boxLength;
            }
        }
    }

    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeTotalMomentum(int N) {
    double numberOfAtoms = 4*N*N*N;
    vec3 totalVelocity = vec3(0, 0, 0);
    for (int i = 0; i < numberOfAtoms; i++){
        Atom *atom = m_atoms[i];
        totalVelocity += atom->velocity;
    }

    for (int i = 0; i < numberOfAtoms; i++){
        Atom *atom = m_atoms[i];
        atom->velocity -= totalVelocity/m_atoms.size();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature, int T_index) {

    double b = UnitConverter::lengthFromAngstroms(latticeConstant);
    m_initialTindex = T_index;

    double N = numberOfUnitCellsEachDimension;
    //Random::randomSeed();
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){

                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));


                atom1->position.set(i*b, j*b, k*b);
                atom2->position.set((0.5+i)*b, (0.5+j)*b, k*b);
                atom3->position.set(i*b, (0.5+j)*b, (0.5+k)*b);
                atom4->position.set((0.5+i)*b, j*b, (0.5+k)*b);

                atom1->initial_position = vec3(i*b, j*b, k*b);
                atom2->initial_position = vec3((0.5+i)*b, (0.5+j)*b, k*b);
                atom3->initial_position = vec3(i*b, (0.5+j)*b, (0.5+k)*b);
                atom4->initial_position = vec3((0.5+i)*b, j*b, (0.5+k)*b);

                atom1->resetVelocityMaxwellian(temperature);
                atom2->resetVelocityMaxwellian(temperature);
                atom3->resetVelocityMaxwellian(temperature);
                atom4->resetVelocityMaxwellian(temperature);
                /*
                double vx1 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vy1 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vz1 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vx2 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vy2 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vz2 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vx3 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vy3 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vz3 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vx4 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vy4 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));
                double vz4 = Random::nextGaussian(0, sqrt(temperature/atom1->mass()));

                atom1->velocity.set(vx1, vy1, vz1);
                atom2->velocity.set(vx2, vy2, vz2);
                atom3->velocity.set(vx3, vy3, vz3);
                atom4->velocity.set(vx4, vy4, vz4);
                */
                m_atoms.push_back(atom1);
                m_atoms.push_back(atom2);
                m_atoms.push_back(atom3);
                m_atoms.push_back(atom4);
            }}}

    setSystemSize(vec3(N*b, N*b, N*b)); // Remember to set the correct system size!
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }

    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
    }



void System::step(double dt, int N) {
    m_integrator.integrate(*this, dt, N);
    m_steps++;
    m_time += dt;
}
