#pragma once 
#include <array>
#include <vector> 
#include <cmath>
#include <utility>
//#include "subcell.hpp"
//#include "SAN_elecphysio.hpp"
//SAN_elecphysio Cell ;
class Mito {
public:
    static constexpr double Km_uni_scale = 4 ; //scaling k constant for MCU (tuned at 6) 
    static constexpr double vNCX_scale = 4 ; //scaling for NCX max transport rate, trying 2/5/10/20/30
    static constexpr double MCU_block_scale = 1 ; //Set equal to % of MCU activity you want (Ru360) - 0.1 = 90% block
    static constexpr double mNCX_block_scale = 1; //Set equl to % of mNCX activity you want (CGP-13757) - 0.1 = 90% block 1
    //Uniporter - Song supplement 3.1 
    static constexpr double Bm_mito = 0.01 ; //Assume MPTP closed - ORIGINAL
    //static constexpr double Bm_mito = 0.03 ; //Try differentvalue 
    //static constexpr double Bm_mito = 0.05 ; //Tuning - Best so far 
    //static constexpr double Bm_mito = 0.08 ; //Tuning
    static constexpr double RTzF = (8.315 * 308) / (96.5 * 2); //z for Ca of +2  
    static constexpr double po_mito = 0.9 ; //MCU open probability 
    static constexpr int N_mcu = 200; //Number of MCUs per mito 
    static constexpr double gMCU_max = 8.1  * MCU_block_scale; //* gMCU_scale //pS Original Value
    static constexpr double Km_uni = 19000 / Km_uni_scale; //micromolar
    static constexpr double zFvmyo = 2.0 * 96.5 * 18 ; //z of 2 for Ca, 18 microl for vmyo

    //Mito NCX - Song supplement 3.2 
    static constexpr double vNCX_max = 0.0035 * vNCX_scale * mNCX_block_scale; //microM * ms^-1, original
    static constexpr double mito_v0 = 91 ; //mV
    static constexpr double kNa_mNCX = 9.4 ; //mM  
    static constexpr double n_mNCX = 3 ; //NCX cooperativity for Nai
    static constexpr double kCa_mNCX = 0.375 ; //micromolar 
    static constexpr double bfRT = (0.5 * 96.5) / (8.316 * 308); //constant b of 0.5 
    static constexpr double nai_mito =  10.0 ; //clamped at 10 
    //double nai_mito =  Cell.y[34]; //clamped na standin 

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
    //static constexpr double DATP = 0.00001 ; //Test
    static constexpr double DATP = 0.01; //Value that yield most heterogeneity (~13 uM difference across CRU grid at 24:1 CRU:mito ratio)
    //static constexpr double DATP = 0.05; 
    //static cons//Mito formulas
    //scaling factors 
    //static constexpr double DATP = 0.25 ; ////diffusion coefficient for ATP (Song)
    static constexpr double l_L_atp = 1.84 ; //longitudinal (x) length constant, um (Song), delta_x
    static constexpr double l_T_atp = 0.9 ; //transverse (y,z) length constant, um (Song), delta_y

    //Compute MCU Flux 
    std::pair<double, double>update_MCU(double ca_space, double mito_psi, double mito_ca){
    double ECa_m = RTzF * log(ca_space/mito_ca) ;
    //double iMCU = (gMCU_max / (1 + (Km_uni/ca_space))) *  (mito_psi - ECa_m);
    double iMCU = (gMCU_max / (1 + (Km_uni/ca_space))) *  (mito_psi + ECa_m); //matching Song's code, sign flip for psi + ECa_m intead of psi - ECa_m
    double J_uni = po_mito * N_mcu * (iMCU/zFvmyo); // Original
    return {iMCU, J_uni} ;
    }

    //Compute mito NCX Flux 
    double update_NCX_mito(double cai, double mito_psi, double mito_ca){
    double numerator = (exp(bfRT * (mito_psi - mito_v0)) * (mito_ca / cai)) ;
    double denominator = pow(1 + kNa_mNCX/nai_mito,n_mNCX) * (1 + (kCa_mNCX / mito_ca)) ;
    double jNCX_m = vNCX_max * (numerator/denominator) ; 
    return jNCX_m ;
    };

    //Compute ATP Production and Consumption 
    double update_ATP_consumption(double ATP, double ADP){
    double psp = ATP/ADP; 
    double VATP_consum = k_ATPconsum * (psp / (psp + k_dpspconsum));
    return VATP_consum;
    }

    double update_ATP_production(double mito_psi, double ATP, double ADP){
    double psp = ATP/ADP; 
    double fdmito = mito_psi / (mito_psi + kd_deltam); 
    double fADP = 1 / (1 + kd_psn * psp);
    double fATP = psp / (psp + kd_psp); 
    double gdmito = 0.3 * (1 - fdmito);
    double VATPase = k_ATPase * (fdmito*fADP - gdmito * fATP) ;
    return VATPase;
    }
};
