
//                dimensionedScalar UxRESV("UxRESV",dimensionSet(0,0,0,0,0,0,0),1e-18);
//                dimensionedScalar UyRESV("UyRESV",dimensionSet(0,0,0,0,0,0,0),1e-18);
//                dimensionedScalar UzRESV("UzRESV",dimensionSet(0,0,0,0,0,0,0),1e-18);
//                UxRESV = max(mag(UxRES - U.component(0))/max(max(mag(U.component(0))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15)));
//                UyRESV = max(mag(UyRES - U.component(1))/max(max(mag(U.component(1))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15)));
//                UzRESV = max(mag(UzRES - U.component(2))/max(max(mag(U.component(2))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15)));
//                dimensionedScalar PRESV("URESV",dimensionSet(0,0,0,0,0,0,0),1e-18);
//                PRESV = max(mag(PRES - p)/max(mag(p),dimensionedScalar("1e-10",dimensionSet(1,-1,-2,0,0,0,0),1e-10)));
//                dimensionedScalar hsRESV("URESV",dimensionSet(0,0,0,0,0,0,0),1e-18);
//                hsRESV = max(mag(hsRES - thermo.T())/max(mag(thermo.T()),dimensionedScalar("1e-10",dimensionSet(0,0,0,1,0,0,0),1e-10)));
//                if(uxres <= max(mag(UxRES - U.component(0))/max(max(mag(U.component(0))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15))).value())
//                {
//                    uxres=max(mag(UxRES - U.component(0))/max(max(mag(U.component(0))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15))).value();
//                }
//                if(uyres <= max(mag(UyRES - U.component(1))/max(max(mag(U.component(1))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15))).value())
//                {
//                    uyres=max(mag(UyRES - U.component(1))/max(max(mag(U.component(1))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15))).value();
//                }
//                if(uzres <= max(mag(UzRES - U.component(2))/max(max(mag(U.component(2))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15))).value())
//                {
//                    uzres=max(mag(UzRES - U.component(2))/max(max(mag(U.component(2))),dimensionedScalar("1e-10",dimensionSet(0,1,-1,0,0,0,0),1e-15))).value();
//                }
//                if(pres <= max(mag(PRES - p)/max(mag(p),dimensionedScalar("1e-10",dimensionSet(1,-1,-2,0,0,0,0),1e-10))).value())
//                {
//                    pres=max(mag(PRES - p)/max(mag(p),dimensionedScalar("1e-10",dimensionSet(1,-1,-2,0,0,0,0),1e-10))).value();
//                }
                if(hsres <= max(mag(hsRES - thermo.T())/max(mag(thermo.T()),dimensionedScalar("1e-10",dimensionSet(0,0,0,1,0,0,0),1e-10))).value())
                {
                    hsres=max(mag(hsRES - thermo.T())/max(mag(thermo.T()),dimensionedScalar("1e-10",dimensionSet(0,0,0,1,0,0,0),1e-10))).value();
                }

