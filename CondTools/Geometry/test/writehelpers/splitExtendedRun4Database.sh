#!/bin/sh

conddb -f sqlite_file:myfile.db -c sqlite_file:GeometryFileExtendedRun4D110.db -t XMLFILE_Geometry_TagXX_ExtendedRun4D110_mc -i XMLFILE_Geometry_TagXX_ExtendedRun4D110_mc
conddb -f sqlite_file:myfile.db -c sqlite_file:TKRECO_Geometry.db -t TKRECO_Geometry_TagXX -i TKRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:TKParameters_Geometry.db -t TKParameters_Geometry_TagXX -i TKParameters_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:EBRECO_Geometry.db -t EBRECO_Geometry_TagXX -i EBRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:HCALRECO_Geometry.db -t HCALRECO_Geometry_TagXX -i HCALRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:HCALParameters_Geometry.db -t HCALParameters_Geometry_TagXX -i HCALParameters_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:HGCALRECO_Geometry.db -t HGCALRECO_Geometry_TagXX -i HGCALRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:HGCALParameters_Geometry.db -t HGCALParameters_Geometry_TagXX -i HGCALParameters_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:CTRECO_Geometry.db -t CTRECO_Geometry_TagXX -i CTRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:ZDCRECO_Geometry.db -t ZDCRECO_Geometry_TagXX -i ZDCRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:CASTORRECO_Geometry.db -t CASTORRECO_Geometry_TagXX -i CASTORRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:CSCRECO_Geometry.db -t CSCRECO_Geometry_TagXX -i CSCRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:CSCRECODIGI_Geometry.db -t CSCRECODIGI_Geometry_TagXX -i CSCRECODIGI_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:DTRECO_Geometry.db -t DTRECO_Geometry_TagXX -i DTRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:RPCRECO_Geometry.db -t RPCRECO_Geometry_TagXX -i RPCRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:GEMRECO_Geometry.db -t GEMRECO_Geometry_TagXX -i GEMRECO_Geometry_TagXX
conddb -f sqlite_file:myfile.db -c sqlite_file:ME0RECO_Geometry.db -t ME0RECO_Geometry_TagXX -i ME0RECO_Geometry_TagXX
