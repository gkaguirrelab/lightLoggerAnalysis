%FLIC_2001 Gaze cal values:
startTime = [1, 23, 683]
targetDurSec = 3.267;
onset_delay_s = 0;
observerArgs = {'sphericalAmetropia',-1.25,'spectacleLens',[-1.25,0,0]};

%%FLIC_2002 Gaze cal values:
startTime = [1, 20, 933];
fullFrameSet = [...
        9812
       10228
       10592
       11437
       11813
       12247
       12605
       13065
       13506
       13776
       13912
       14435
       14791
       15168
       15644
       16033
       16470
       16937
       17289
       17686
       18038
       18609
       18974
       19308
       19779
       20225
       20692
       21039
       21510
       21857
       21995
       22395
       22809]; % 13776; 21995 had to be added manually because these targets were skipped. two extra frames from the end were deleted. Target 4 had bad perimeters so we skipped target 4.
% skip target 4
gazeTargets =

         0         0
   15.0000   15.0000
   15.0000  -15.0000
  -15.0000  -15.0000
         0   15.0000
         0  -15.0000
   15.0000         0
  -15.0000         0
    7.5000    7.5000
    7.5000   -7.5000
   -7.5000    7.5000
   -7.5000   -7.5000
         0   10.0000
         0   -7.5000
    7.5000         0
   -7.5000         0
         0         0
   15.0000   15.0000
   15.0000  -15.0000
  -15.0000   15.0000
  -15.0000  -15.0000
         0   15.0000
         0  -15.0000
   15.0000         0
  -15.0000         0
    7.5000    7.5000
    7.5000   -7.5000
   -7.5000    7.5000
   -7.5000   -7.5000
         0   10.0000
         0   -7.5000
    7.5000         0
   -7.5000         0

   observerArgs = {'sphericalAmetropia',-1.25};
p5 =

  Columns 1 through 7

  -28.6484   -7.3094   51.0564   24.3158    0.5042   12.1706    0.9918

  Columns 8 through 11

    0.9927   18.8754   49.3395   40.5355

pMean = [-20.2767
  -11.7980
   56.2349
   16.5126
   -4.3540
   14.9336
    1.0004
    1.0328
   12.9956
   49.0111
   40.2369];

   p34 = [-20.2778  -11.7970   56.2308   16.5185   -4.3553   14.9316    1.0004 1.0328   12.4987   49.9999   40.0000];
    
gazeOffset = [-0.8, 0.3] % [azi, ele]
%NOTE! frame selection was done with default of 0.8 confidence for at least
%8 points. scene geometry was estimated with default of 0.75 confidence for
%each point.


%%FLIC_2003
startTime = [1, 24, 525]; % [minutes, seconds, milliseconds]
firstDotEnd = [1 26 142];
secondDotEnd = [1 39 392];% actually 5th probably because of a long buffering gap
thirdDotEnd =[1 42 667];
targetDurSec = 3.3;
% the gaps are very long for this participant and the confidence was too
% low on many point
confidenceCutoff = 0.6;
fullFrameSet =[...
       10223
       10552
       10899
       11671
       12122
       12553
       12986
       13265
       13758
       14046
       14425
       14882
       15309
       15724
       16081
       16415
       16958
       17271
       17744
       18014
       18442
       18826
       19282
       19683
       20095
       20422
       20808
       21178
       21634
       22088
       22489
       22805
       23201];

gazeTargetsDeg =

         0         0
  -15.0000   15.0000
  -15.0000  -15.0000
   15.0000  -15.0000
         0   15.0000
         0  -15.0000
  -15.0000         0
   15.0000         0
   -7.5000    7.5000
   -7.5000   -7.5000
    7.5000    7.5000
    7.5000   -7.5000
         0    7.5000
         0   -7.5000
   -7.5000         0
    7.5000         0
         0         0
  -15.0000   15.0000
  -15.0000  -15.0000
   15.0000   15.0000
   15.0000  -15.0000
         0   15.0000
         0  -15.0000
  -15.0000         0
   15.0000         0
   -7.5000    7.5000
   -7.5000   -7.5000
    7.5000    7.5000
    7.5000   -7.5000
         0    7.5000
         0   -7.5000
   -7.5000         0
    7.5000         0

observerArgs = {'sphericalAmetropia',-5.75,'spectacleLens',[-4.5,0,0]};


%%FLIC_2004
startTime = [1, 39, 408]; % [minutes, seconds, milliseconds]
firstDotEnd = [1 41 908];
secondDotEnd = [1 45 275];
thirdDotEnd =[1 48 608];
targetDurSec = 3.600;
confidenceCutoff = 0.8;

fullFrameSet =[...

       12010
       12339
       12756
       13226
       14086
       14868
       15218
       15733
       16117
       16589
       16976
       17461
       17882
       18329
       18763
       19155
       20137
       20954
       21370
       21804
       22239
       22621
       23076
       23508
       23946
       24386
       24813
       25156
       25658
       26173];

% had to pause on FLIC_2004 because we need optical prescription

%%FLIC_2005
startTime = [1, 25, 400]; % [minutes, seconds, milliseconds]
firstDotEnd = [1 26 608];
secondDotEnd = [1 30 000]; %actually had to skip because of a buffering section
thirdDotEnd =[1 39 983];

goodIdx =

     1
     2
     4
     5
     6
     7
     8
     9
    10
    11
    12
    13
    14
    15
    16
    17
    18
    19
    20
    21
    22
    23
    24
    25
    26
    27
    28
    29
    30
    31
    32
    33
    34


fullFrameSet = [...
       10349
       10666
       11507
       11929
       12362
       12783
       13153
       13518
       13962
       14319
       14661
       15173
       15494
       15869
       16275
       16656
       17154
       17502
       17951
       18348
       18786
       19056
       19494
       19920
       20391
       20763
       21072
       21545
       21951
       22327
       22741
       23083
       23593]; %skipped target 3 which was mostly during a buffering break
gazeTargetsDeg =

         0         0
  -15.0000   15.0000
   15.0000   15.0000
   15.0000  -15.0000
         0   15.0000
         0  -15.0000
  -15.0000         0
   15.0000         0
   -7.5000    7.5000
   -7.5000   -7.5000
    7.5000    7.5000
    7.5000   -7.5000
         0    7.5000
         0   -7.5000
   -7.5000         0
    7.5000         0
         0         0
  -15.0000   15.0000
  -15.0000  -15.0000
   15.0000   15.0000
   15.0000  -15.0000
         0   15.0000
         0  -15.0000
  -15.0000         0
   15.0000         0
   -7.5000    7.5000
   -7.5000   -7.5000
    7.5000    7.5000
    7.5000   -7.5000
         0    7.5000
         0   -7.5000
   -7.5000         0
    7.5000         0

    observerArgs = {'sphericalAmetropia',-1.00,'spectacleLens',[-1.00,-0.25,0]};
