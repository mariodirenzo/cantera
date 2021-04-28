//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/IonFlow.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/numerics/polyfit.h"

using namespace std;

namespace Cantera
{

IonFlow::IonFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    StFlow(ph, nsp, points),
    m_import_electron_transport(false),
    m_stage(1),
    m_kElectron(npos)
{
    // make a local copy of species charge
    for (size_t k = 0; k < m_nsp; k++) {
        m_speciesCharge.push_back(m_thermo->charge(k));
    }

    // Find indices for charge of species
    for (size_t k = 0; k < m_nsp; k++){
        if (m_speciesCharge[k] != 0){
            m_kCharge.push_back(k);
        } else {
            m_kNeutral.push_back(k);
        }
    }

    // Find the index of electron
    if (m_thermo->speciesIndex("E") != npos ) {
        m_kElectron = m_thermo->speciesIndex("E");
    }

    // no bound for electric potential
    setBounds(c_offset_P, -1.0e20, 1.0e20);

    // Set tighter negative species limit on charged species to avoid
    // instabilities. Tolerance on electrons is even tighter to account for the
    // low "molecular" weight.
    for (size_t k : m_kCharge) {
        setBounds(c_offset_Y + k, -1e-8, 1.0);
    }
    if (m_kElectron != npos) {
        setBounds(c_offset_Y + m_kElectron, -1e-10, 1.0);
    }

    m_refiner->setActive(c_offset_P, true);
    m_mobility.resize(m_nsp*m_points);
    m_do_poisson.resize(m_points,false);
}

void IonFlow::resize(size_t components, size_t points){
    StFlow::resize(components, points);
    m_mobility.resize(m_nsp*m_points);
    m_do_species.resize(m_nsp,true);
    m_do_poisson.resize(m_points,false);
    m_fixedElecPoten.resize(m_points,0.0);
}

void IonFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x,j0,j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobility[j*m_nsp]);
        if (m_import_electron_transport) {
            size_t k = m_kElectron;
            double tlog = log(m_thermo->temperature());
            m_mobility[k+m_nsp*j] = poly5(tlog, m_mobi_e_fix.data());
            m_diff[k+m_nsp*j] = poly5(tlog, m_diff_e_fix.data());
        }
    }
}

void IonFlow::updateDiffFluxes(const double* x, size_t j0, size_t j1)
{
    if (m_stage == 1) {
        frozenIonMethod(x,j0,j1);
    }
    if (m_stage == 2) {
        poissonEqnMethod(x,j0,j1);
    }
}

void IonFlow::frozenIonMethod(const double* x, size_t j0, size_t j1)
{
    StFlow::updateDiffFluxes(x, j0, j1);
//    for (size_t j = j0; j < j1; j++) {
//        // flux for ions
//        // Set flux to zero to prevent some fast charged species (e.g. electron)
//        // to run away
//        for (size_t k : m_kCharge) {
//            m_flux(k,j) = 0;
//        }
//    }
}

void IonFlow::poissonEqnMethod(const double* x, size_t j0, size_t j1)
{
    // Add diffusion fluxes
    StFlow::updateDiffFluxes(x, j0, j1);

    // Add ambipolar diffusion
    for (size_t j = j0; j < j1; j++) {
        const double E_ambi = E(x,j);
        const double rho = 0.5*(density(j) + density(j+1));
        double sum = 0.0;
        for (size_t k : m_kCharge) {
            const double Vdrift = m_speciesCharge[k] * m_mobility[k+m_nsp*j] * E_ambi;
            // Upwind the mass fraction reconstruction based on the sign of drift velocity
            const double Flux = rho * Vdrift * ((Vdrift > 0.0) ? Y(x,k,j) : Y(x,k,j+1));
            m_flux(k,j) += Flux;
            sum -= Flux;
        }
        // correction flux to insure that \sum_k Y_k V_k = 0.
        for (size_t k = 0; k < m_nsp; k++) {
            m_flux(k,j) += sum*Y(x,k,j);
        }
    }
}

void IonFlow::setSolvingStage(const size_t stage)
{
    if (stage == 1 || stage == 2) {
        m_stage = stage;
    } else {
        throw CanteraError("IonFlow::setSolvingStage",
                    "solution stage must be set to: "
                    "1) frozenIonMethod, "
                    "2) poissonEqnMethod");
    }
}

void IonFlow::evalResidual(double* x, double* rsd, int* diag,
                           double rdt, size_t jmin, size_t jmax)
{
    StFlow::evalResidual(x, rsd, diag, rdt, jmin, jmax);
    for (size_t j = jmin; j <= jmax; j++) {
        if (j == 0) {
            // Boundary conditions will take care of this point later
            rsd[index(c_offset_P, j)] = phi(x,j);
            diag[index(c_offset_P, j)] = 0;
        } else if (j == m_points - 1) {
            // Boundary conditions will take care of this point later
            rsd[index(c_offset_P, j)] = phi(x,j);
            diag[index(c_offset_P, j)] = 0;
        } else {
            //-----------------------------------------------
            //    Energy equation
            //    Add Jule heating
            //
            //    Na*sum(q_k*J_k/W_k)*E
            //-----------------------------------------------
            if (m_do_energy[j]) {
                double JuleHeating = 0.0;
                for (size_t k : m_kCharge) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    JuleHeating += m_speciesCharge[k] * flxk / m_wt[k];
                }
                JuleHeating *= ElectronCharge * Avogadro * 0.5*(E(x,j) + E(x,j-1));
                rsd[index(c_offset_T, j)] += JuleHeating/(m_rho[j]*m_cp[j]);
            }
            //-----------------------------------------------
            //    Poisson's equation
            //
            //    dE/dz = e/eps_0 * sum(q_k*n_k)
            //
            //    E = -dPhi/dz
            //-----------------------------------------------
            if (m_stage == 2) {
               rsd[index(c_offset_P, j)] = dEdz(x,j) - rho_e(x,j)/epsilon_0;
            } else {
               rsd[index(c_offset_P, j)] = dEdz(x,j);
            }
            diag[index(c_offset_P, j)] = 0;
        }
    }
}

void IonFlow::solvePoissonEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = true;
        }
    } else {
        if (!m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    m_refiner->setActive(c_offset_P, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixElectricPotential(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = false;
        }
    } else {
        if (m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_P, false);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::setElectronTransport(vector_fp& tfix, vector_fp& diff_e,
                                   vector_fp& mobi_e)
{
    m_import_electron_transport = true;
    size_t degree = 5;
    size_t n = tfix.size();
    vector_fp tlog;
    for (size_t i = 0; i < n; i++) {
        tlog.push_back(log(tfix[i]));
    }
    vector_fp w(n, -1.0);
    m_diff_e_fix.resize(degree + 1);
    m_mobi_e_fix.resize(degree + 1);
    polyfit(n, degree, tlog.data(), diff_e.data(), w.data(), m_diff_e_fix.data());
    polyfit(n, degree, tlog.data(), mobi_e.data(), w.data(), m_mobi_e_fix.data());
}

void IonFlow::_finalize(const double* x)
{
    StFlow::_finalize(x);

    //bool p = m_do_poisson[0];
    //for (size_t j = 0; j < m_points; j++) {
    //    if (!p) {
    //        m_fixedElecPoten[j] = phi(x, j);
    //    }
    //}
    //if (p) {
    //    solvePoissonEqn();
    //}
}

}
