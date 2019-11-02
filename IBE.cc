/*

This software is a computer program which purpose is to provide 
a proof-of-concept implementation of PKC 2020 submission
"A New Trapdoor over Module-NTRU Lattices and its Application to ID-based Encryption".

This software is based on Thomas Prest's original implementation
governed by the CeCILL license under French law and abiding by the rules of distribution of free software.
The original implementation can be found in "https://github.com/tprest/Lattice-IBE/"

*/



#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>


#include "params.h"
#include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "Scheme.h"


using namespace std;
using namespace NTL;




//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================


int main()
{   
    ZZ_pX phiq, MPK;
    unsigned int i;
    float diff;
    clock_t t1, t2;
    const ZZX phi = Cyclo();

    srand(rdtsc()); // initialisation of rand

    cout << "N = " << N0 << endl;
    cout << "q = " << q0 << endl;

    ZZ_p::init(q1);
    zz_p::init(q0);

    phiq = conv<ZZ_pX>(phi);
    ZZ_pXModulus PHI(phiq);

    // NTRU
    // cout << "\n=======================================================================\n";
    // cout << "This program is a proof-of concept for efficient IBE over NTRU lattices.\n";
    // cout << "It generates a NTRU lattice of dimension 2N and associated modulus q,\n";
    // cout << "and perform benches and tests, for user key extraction and encryption/decryption.";
    // cout << "\n=======================================================================\n\n";

    // ZZX MSK[4];
    // MSK_Data * MSKD = new MSK_Data;
    // MPK_Data * MPKD = new MPK_Data;



    // cout << "\n===================================================================\n KEY GENERATION";
    // cout << "\n===================================================================\n";
    // t1 = clock();
    // //NTRU
    // for(i = 0; i < 1; i++)
    // {
    //     Keygen(MPK, MSK);
    // }

    // CompleteMSK(MSKD, MSK);
    // CompleteMPK(MPKD, MPK);

    // cout << "GS norm of MSK = max(" << MSKD->GS_Norms[0] << ", " << MSKD->GS_Norms[N0] << ")" <<endl;


    // t2 = clock();
    // diff = ((float)t2 - (float)t1)/1000000.0F;
    // cout << "It took " << diff << " seconds to generate the Master Secret Key" << endl;



    // //==============================================================================
    // //Key extraction bench and encryption/decryption bench
    // //==============================================================================
    // const unsigned int nb_extrb = 100;
    // const unsigned int nb_crypb = 1000;

    // cout << "\n===================================================================\n RUNNING EXTRACTION BENCH FOR ";
    // cout << nb_extrb << " DIFFERENT IDENTITIES\n===================================================================\n";
    // Extract_Bench(nb_extrb, MSKD);

    // cout << "\n===================================================================\n RUNNING ENCRYPTION BENCH FOR ";
    // cout << nb_crypb << " DIFFERENT MESSAGES\n===================================================================\n";
    // Encrypt_Bench(nb_crypb, MPKD, MSKD);


    // ///==============================================================================
    // //Key extraction test and encryption/decryption test
    // //==============================================================================
    // const unsigned int nb_extrt = 100;
    // const unsigned int nb_crypt = 100;

    // cout << "\n===================================================================\n CHECKING EXTRACTION VALIDITY FOR ";
    // cout << nb_extrt << " DIFFERENT IDENTITIES\n===================================================================\n";
    // Extract_Test(nb_extrt, MSKD);

    // cout << "\n===================================================================\n CHECKING ENCRYPTION VALIDITY FOR ";
    // cout << nb_extrt << " DIFFERENT MESSAGES\n===================================================================\n";
    // Encrypt_Test(nb_crypt, MPKD, MSKD);

    // free(MSKD);
    // free(MPKD);
    // return 0;

//====================================================================================================
//====================================================================================================
//====================================================================================================
//====================================================================================================    
//GNTRU
    cout << "\n=======================================================================\n";
    cout << "This program is a proof-of concept for efficient IBE over GNTRU lattices.\n";
    cout << "It generates a GNTRU lattice of dimension 3N and associated modulus q,\n";
    cout << "and perform benches and tests, for user key extraction and encryption/decryption.";
    cout << "\n=======================================================================\n\n";
    cout << "\n===================================================================\n KEY GENERATION";
    cout << "\n===================================================================\n";

    ZZX GNTRU_MSK[9];
    ZZ_pX GNTRU_MPK[2];
    GNTRU_MSK_Data * GNTRU_MSKD = new GNTRU_MSK_Data;
    GNTRU_MPK_Data * GNTRU_MPKD = new GNTRU_MPK_Data;

    t1 = clock();
    for(i = 0; i < 1; i++)
    {
        GNTRU_Keygen(GNTRU_MPK, GNTRU_MSK);
    }

    GNTRU_CompleteMSK(GNTRU_MSKD, GNTRU_MSK);
    GNTRU_CompleteMPK(GNTRU_MPKD, GNTRU_MPK);

    cout << "GS norm of MSK = max(" << GNTRU_MSKD->GS_Norms[0] << ", " << GNTRU_MSKD->GS_Norms[N0] << ", " << GNTRU_MSKD->GS_Norms[2 * N0] << ")" <<endl;

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "It took " << diff << " seconds to generate the Master Secret Key" << endl;

    //==============================================================================
    //Key extraction bench and encryption/decryption bench
    //==============================================================================
    const unsigned int GNTRU_nb_extrb = 100;
    const unsigned int GNTRU_nb_crypb = 1000;

    cout << "\n===================================================================\n RUNNING EXTRACTION BENCH FOR ";
    cout << GNTRU_nb_extrb << " DIFFERENT IDENTITIES\n===================================================================\n";
    GNTRU_Extract_Bench(GNTRU_nb_extrb, GNTRU_MSKD);

    cout << "\n===================================================================\n RUNNING ENCRYPTION BENCH FOR ";
    cout << GNTRU_nb_crypb << " DIFFERENT MESSAGES\n===================================================================\n";
    GNTRU_Encrypt_Bench(GNTRU_nb_crypb, GNTRU_MPKD, GNTRU_MSKD);


    ///==============================================================================
    //Key extraction test and encryption/decryption test
    //==============================================================================
    const unsigned int GNTRU_nb_extrt = 100;
    const unsigned int GNTRU_nb_crypt = 100;

    cout << "\n===================================================================\n CHECKING EXTRACTION VALIDITY FOR ";
    cout << GNTRU_nb_extrt << " DIFFERENT IDENTITIES\n===================================================================\n";
    GNTRU_Extract_Test(GNTRU_nb_extrt, GNTRU_MPKD, GNTRU_MSKD);

    cout << "\n===================================================================\n CHECKING ENCRYPTION VALIDITY FOR ";
    cout << GNTRU_nb_extrt << " DIFFERENT MESSAGES\n===================================================================\n";
    GNTRU_Encrypt_Test(GNTRU_nb_crypt, GNTRU_MPKD, GNTRU_MSKD);
  

    free(GNTRU_MSKD);
    free(GNTRU_MPKD);
    return 0;
}
