// RELATIVISTIC KINEMATICS (ROOT program)
// Program for relativistic kinematic calculations (elastic and inelastic)
// Programmed by R.E. Philips and S.T. Thornton
// Oak Ridge national lab, july 1967.
// Reload by Andreyanov A.

///////////////////////////////////////////////////////////////////////
// Reaction 1 + 2 --> 3 + 4

// c = 1

// INPUT DATA
// 1) m1, m2, m3, m4 -- particles masses of rest (Mev); 
// 2) E1_k -- kinetic energy of particle 1 in laboratory system;
// 3) theta -- laboratory scattering angle of particle 3.

// OUTPUT DATA
// 1) E1_cm; E2_cm; E3_cm; E4_cm; 
// 2) Theta3_cm; Theta4_cm; Theta4;
// 3) E3; E4.

#include <TMath.h>
using namespace std;

void Relativ_kinematic_calculations(){
    double m1, m2, m3, m4, Q, E1_k, theta, EN, Q_gs, T, E1, E_th, ET, p1,\
    T_A, sin_th3, cos_th3, check, E_cmT, E1_cm, E2_cm, E3_cm,  E4_cm, \
    p_cm, alpha, E3, E4, p3, p4, beta, gamma, cos_th3_cm, sin_th3_cm, \
    th3_cm, th4_cm, sin_th4, cos_th4, theta4, E1_k_min;

    // Input data (Mev and degrees)
    m1 = 1875.61;
    m2 = 1875.61;
    m3 = 2808.92;
    m4 = 939.565;
    E1_k = 100;
    theta = 5;

    Q_gs = m1 + m2 - m3 - m4; // Q-value groundstate
    E1_k_min = (-Q_gs*(m1+m2+m3+m4)/(2*m2)); // Threshold energy
    // Check it is enough energy for this reaction
    if (E1_k <= E1_k_min){
        cout << "Not enough kinetic energy for this reaction" << "\n";
        return;
    }
    E1 = E1_k + m1; // Full energy of particle 1
    E_th = E1_k_min + m1; // Full threshold energy
    ET = E1 + m2; // Total energy (ET = E1+E2 = E3+E4)
    p1 = sqrt(E1*E1 - m1*m1); //  Momentum of particle 1 in lab
    
    //Change angles to radians
    theta = theta*0.17453292;
    sin_th3 = sin(theta);
    cos_th3 = cos(theta);

    // Check it is enough energy for this angle
    check = pow(m2*E1 + (m1*m1 + m2*m2 - m3*m3 -m4*m4)/2,2) \
    - pow(m3*m4,2) - pow(p1*m3*sin_th3,2);
    //if (check < 0){
    //    cout << "Not enough energy for this angle " << "\n";
    //    return;
    //}

    // Calculate CENTER OF MASS energies
    E_cmT = sqrt(m1*m1 + m2*m2 + 2*m2*E1);
    E1_cm = (m1*m1 + m2*E1)/E_cmT;
    E2_cm = (m2*m2 + m2*E1)/E_cmT;
    E3_cm = (E_cmT*E_cmT + m3*m3 - m4*m4)/(2*E_cmT);
    E4_cm = (E_cmT*E_cmT + m4*m4 - m3*m3)/(2*E_cmT);

    p_cm = sqrt(E4_cm*E4_cm - m4*m4); // Center of mass momentum 
                                      // of outgoing particles 

    // Alpha determines if there two or one solutions for E3
    alpha = (p1*(1+(m3*m3 - m4*m4)/(E_cmT*E_cmT)))/\
    (ET*sqrt((1-pow((m3 + m4)/E_cmT,2))*(1-pow((m3-m4)/E_cmT,2))));
    if (alpha<=0){
        cout << "One root for E3" << "\n";
        const int n=1;
        double E3[n];
        E3[0] = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2)\
        + p1*cos_th3*sqrt(check))/(ET*ET - pow(p1*cos_th3,2));
    } else{
        cout << "Two roots for E3" << "\n";
        const int n=2;
        double E3[n];
        E3[0] = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2)\
        + p1*cos_th3*sqrt(check))/(ET*ET - pow(p1*cos_th3,2));
        E3[1] = (ET*(m2*E1 + (m1*m1 +m2*m2 + m3*m3 - m4*m4)/2)\
        - p1*cos_th3*sqrt(check))/(ET*ET - pow(p1*cos_th3,2));
    }

    double E4[n], p3[n], p4[n], sin_th3_cm[n], cos_th3_cm[n], th3_cm[n], th4_cm[n],\
    sin_th4[n], cos_th4[n], theta4[n], dsigma_cm__dsigma_lab[n], 
    beta = p1/ET;
    gamma = ET/E_cmT;
    // Make loop for finding roots: one step for alpha < 0 and 
    // two steps for alpha > 0
    for (int i=0; i<n; i++){
        E4[i] = ET - E3[i];
        p3[i] = sqrt(E3[i]*E3[i] - m3*m3);
        p4[i] = sqrt(E4[i]*E4[i] - m4*m4);
        // Calculate center of mass angles of particles 3 and 4
        sin_th3_cm[i] = p3[i]*sin_th3/p_cm;
        cos_th3_cm[i] = gamma*(p3[i]*cos_th3 - beta*E3[i])/p_cm;
        th3_cm[i] = TMath::ATan(sin_th3_cm[i]/cos_th3_cm[i]); // Theta3_cm
        if (cos_th3_cm[i] < 0){
            th3_cm[i] = 3.14159265 + th3_cm[i];
        }
        th4_cm[i] = 3.14159265 - th3_cm[i]; // In CM system theta3_cm + theta4_cm =180
        // Calculate lab scattering angle of particle 4
        sin_th4[i] = p3[i]*sin_th3/p4[i];
        cos_th4[i] = (p1 - p3[i]*cos_th3)/p4[i];
        theta4[i] = TMath::ATan(sin_th4[i]/cos_th4[i]);
        if (cos_th4[i] < 0){
            theta4[i] = 3.14159265 + theta4[i];
        }

        // Ratio differential cross section in CM to LAB
        dsigma_cm__dsigma_lab[i] = ET/E_cmT*(1 + (1 -2*i)*p1/ET*(ET*ET - p1*p1 + m3*m3 -\
        m4*m4)/sqrt(pow(ET*ET - p1*p1 + m3*m3 - m4*m4, 2) - 4*E_cmT*E_cmT*m3*m3)*\
        cos_th3_cm[i])*pow(sin_th3,3)/pow(sin_th3_cm[i],3);

        // Translate angles to degrees
        th4_cm[i] = th4_cm/0.17453292;
        th3_cm[i] = th3_cm/0.17453292;
        theta4[i] = theta4/0.17453292;

        // Change energies to kinetic (MeV)
        E3[i] = E3[i] - m3;
        E4[i] = E4[i] - m4;
    }

    // Change energies to kinetic (MeV)
    E1_cm = E1_cm - m1;
    E2_cm = E2_cm - m2;
    E3_cm = E3_cm - m3;
    E4_cm = E4_cm - m4;
    theta = theta/0.17453292;

    // PRINT OUTPUT
    cout << "OUTPUT \nINPUT DATA" << "\n";
    printf("m1 %f; m2 %f; m3 %f; m4 %f; E_kin %f; theta %f \n", m1, m2, m3, m4, E1_k, theta);
    cout << "OUTPUT DATA" << "\n";
    printf("E1_cm %f\n", E1_cm);
    printf("E2_cm %f\n", E2_cm);
    printf("E3_cm %f\n", E3_cm);
    printf("E4_cm %f\n", E4_cm);

    for(int i2=0; i2<n; i2++){
        
        printf("SOLUTION %f\n", i2+1);
        printf("Theta3_cm %f\n", th3_cm[i2]);
        printf("Theta4_cm %f\n", th4_cm[i2]);
        printf("Theta4 %f\n", theta4[i2]);

        printf("E3 %f\n", E3[i2]);
        printf("E4 %f\n", E4[i2]);
    }
   
}
