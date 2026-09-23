#pragma once 
#include <array>
#include <vector> 
#include <cmath>
#include <utility>

//Mito formulas
//Uniporter - Song supplement 3.1 
static constexpr double Bm_mito = 0.01 ; //Assume MPTP closed
static constexpr double RTzF = (8.315 * 308) / (96.5 * 2); //z for Ca of +2  
static constexpr double po_mito = 0.9 ; //MCU open probability 
static constexpr int N_mcu = 200 ; //Number of MCUs per mito 
static constexpr double gMCU_max = 8.1 ; //pS
static constexpr double Km_uni = 19000 ; //micromolar
static constexpr double zFvmyo = 2.0 * 96.5 * 18 ; //z of 2 for Ca, 18 microl for vmyo

//Mito NCX - Song supplement 3.2 
static constexpr double vNCX_max = 0.0035 ; //microM * ms^-1
static constexpr double mito_v0 = 91 ; //mV
static constexpr double kNa_mNCX = 9.4 ; //mM  
static constexpr double n_mNCX = 3 ; //NCX cooperativity for Nai
static constexpr double kCa_mNCX = 0.375 ; //micromolar 
static constexpr double bfRT = (0.5 * 96.5) / (8.316 * 308); //constant b of 0.5 
static constexpr double nai_mito = 10.0 ; //clamped na standin 

//Mito Membrane Potential (psi) - Song Supplement 3.6
static constexpr double V_mitos = 3.5 ; //mV*ms*-1
static constexpr double C_mito = 1.812 ; //microMolar * mv-1
static constexpr double k_mitou = 0.0192 ; //ms-1
static constexpr double z_Ca = 2.0 ; //Valence for Ca 

//ATP constants - Song supplement 4.1 
static constexpr double TAN = 7000 ; //Total nucleotide content
static constexpr double ADP_buffer_rate = 0.025; //Proportion of free cyto ADP 
static constexpr double k_ATPconsum = 0.01 ; 
static constexpr double k_dpspconsum = 200 ; 
static constexpr double k_ATPase = 0.16 ; 
static constexpr double kd_deltam = 150;
static constexpr double kd_psn = 0.02 ; 
static constexpr double kd_psp = 0.5 ; 
static constexpr double DATP = 0.25 ; //diffusion coefficient
static constexpr double dx = 1.84 ; //distance of diffusion
static constexpr double tau_ATP = dx * dx / DATP ; //Time constant for ATP diffusion function

//Compute MCU Flux 
inline std::pair<double, double>update_MCU(double ca_space, double mito_psi, double mito_ca){
double ECa_m = RTzF * log(ca_space/mito_ca) ;
double iMCU = (gMCU_max / (1 + (Km_uni/ca_space))) *  (mito_psi - ECa_m);
double J_uni = po_mito * N_mcu * (iMCU/zFvmyo) ;
return {iMCU, J_uni} ;
}

//Compute mito NCX Flux 
inline double update_NCX_mito(double cai, double mito_psi, double mito_ca){
double numerator = (exp(bfRT * (mito_psi - mito_v0)) * (mito_ca / cai)) ;
double denominator = pow(1 + kNa_mNCX/nai_mito,n_mNCX) * (1 + (kCa_mNCX / mito_ca)) ;
double jNCX_m = vNCX_max * (numerator/denominator) ; 
return jNCX_m ;
};

//Compute ATP Production and Consumption 
inline std::pair<double, double> update_ATP_rates(double mito_psi, double ATP, double ADP){
double psp = ATP/ADP; 
double fdmito = mito_psi / (mito_psi + kd_deltam); 
double fADP = 1 / (1 + kd_psn * psp);
double fATP = psp / (psp + kd_psp); 
double gdmito = 0.3 * (1 - fdmito);
double VATP_consum = k_ATPconsum * (psp / (psp + k_dpspconsum));
double VATPase = k_ATPase * (fdmito*fADP - gdmito * fATP) ;
return {VATP_consum, VATPase};
}
