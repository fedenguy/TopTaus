void teenyWeenyCalcNewFakes(){
  
  // MC exp
  double 
    mcExp(2690.7),
    mcErr(77.3);
  
  
  double
    // Fractions
    a1(0.62),
    b1(0.38),

    a2(0.88),
    b2(0.12),
    
    ai(0.835),
    bi(0.165),

    // Errors
    ea1(0.10*a1),
    ea2(0.10*a2),
    eb1(0.10*b1),
    eb2(0.10*b2),
    eai(0.10*ai),
    ebi(0.10*bi);

  double beta  = ( ai*b1 - bi*a1  ) / ( a2*b1 - b2*a1 );
  double alpha = (ai - beta*a2 )/ a1;
  

  double
    i1(3197.9), 
    i2(2682.0); 

  double
    oldEst( (i1+i2)/2. ),
    oldErr( fabs((i1-i2)/2.) ),
    
    newEst( alpha*i1 + beta*i2 ),
    newErr( fabs(alpha*i1 - beta*i2) );

  double
    newAlpha_a(  (ai - a2)/(a1 - a2) ),
    newAlpha_b(  (bi - b2)/(b1 - b2) );

  
  double 
    N_a( ai - a2 ),
    D_a( a1 - a2 ),
    N_b( bi - b2 ),
    D_b( b1 - b2),
    
    eN_a( eai + ea2 ),
    eD_a( ea1 + ea2 ),
    eN_b( ebi + eb2 ),
    eD_b( eb1 + eb2);
  
  double errNewAlpha_a( sqrt(   eN_a*D_a*eN_a*D_a   + N_a*eD_a*N_a*eD_a  )/D_a*D_a            );
  double errNewAlpha_b( sqrt(   eN_b*D_b*eN_b*D_b   + N_b*eD_b*N_b*eD_b  )/D_b*D_b            );

  double 
    recEst_a( newAlpha_a*i1 + (1-newAlpha_a)*i2  ),
    recErr_a( errNewAlpha_a*i1 + errNewAlpha_a*i2   ),
    
    recEst_b( newAlpha_b*i1 + (1-newAlpha_b)*i2   ),
    recErr_b( errNewAlpha_b*i1 + errNewAlpha_b*i2   );

 
  cout << "MC expected : " << mcExp  << " +/- " << mcErr << endl;
  cout << "---------------------" << endl;
  cout << "Old method" << endl;
  cout << "Old W+j  value: " << i1 << endl;
  cout << "Old QCDm value: " << i2 << endl;
  cout << "Formula for old estimate: " << "0.5*" << i1 << " + " << "0.5*" << i2 << endl;   
  cout << "Old estimate: " << oldEst << " +/- " << oldErr << " (" << 100*oldErr/oldEst << "%) " << endl;
  cout << "---------------------" << endl;
  cout << "a1: " << a1 << endl;
  cout << "b1: " << b1 << endl;
  cout << "a2: " << a2 << endl;
  cout << "b2: " << b2 << endl;
  cout << "ai: " << ai << endl;
  cout << "bi; " << bi << endl;
  cout << "---------------------" << endl;
  cout << "Method 1: alpha*ai + beta*bi " << endl;
  cout << "alpha: " << alpha << endl;
  cout << "beta : " << beta << endl;
  cout << "New W+j  value: " << alpha*i1 << endl;
  cout << "New QCDm value: " << beta*i2  << endl;
  cout << "Formula for new estimate: " << alpha << "*" << i1 << " + " << beta <<  "*" << i2 << endl;
  cout << "New estimate: " << newEst << " +/- " << newErr << " (" << 100*newErr/newEst << "%) " << endl;  
  cout << "---------------------" << endl;
  cout << "Method 2 alpha*ai + (1-alpha)*bi" << endl;
  cout << "alpha (through as): " << newAlpha_a << " +/- " << errNewAlpha_a << " (" << 100*errNewAlpha_a/newAlpha_a << "%)" << endl;
  cout << "alpha (through bx): " << newAlpha_b << " +/- " << errNewAlpha_b << " (" << 100*errNewAlpha_b/newAlpha_b << "%)" << endl;
  cout << "Formula for new estimate: " << newAlpha_a << "*" << i1 << " + " << 1-newAlpha_a <<  "*" << i2 << endl;
  cout << "New estimate (through as): " << recEst_a << " +/- " << recErr_a << " (" << 100*recErr_a/recEst_a << "%) " << endl;
  cout << "New estimate (through bs): " << recEst_b << " +/- " << recErr_b << " (" << 100*recErr_b/recEst_b << "%) " << endl;
  cout << "---------------------" << endl;
  cout << "----------------" << endl;




}
