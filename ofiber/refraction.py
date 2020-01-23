# pylint: disable=invalid-name
# pylint: disable=line-too-long
"""
Useful routines for index of refraction based on Sellmeier coefficients

Included are Sellmeier coefficients for a bunch of common glasses used
in optical fibers.  Functions are provided to generate the index of refraction,
as well as its first and second derivatives.

    Use like this
        import ofiber as of
        lambdas = np.linspace(1000,1700,50)*1e-9 # [m]
        i = of.find_glass("SiO2")
        glass = of.glass(i)
        n = of.n(glass,lambdas)
        plt.plot(lambdas*1e9, n)

Todo:
    * add fancy names
    * throw an exception when name not found?
    * add proper citations
    * convert to objects?
    * move commented stuff at end to documentation

Scott Prahl
Aug 2018
"""

import numpy as np


__all__ = ('ALL_GLASS_NAMES',
           'd2n',
           'dn',
           'doped_glass',
           'doped_glass_name',
           'glass',
           'glass_name',
           'find_glass',
           'n',
           'n_group',
           'n_air')

_GLASS = [
    # format is [c1, c2, c3, b1, b2, b3]
    # where b1, b2, b3 are unitless and c1, c2, c3 have units of [microns**2]
    [4.67914826e-3, 1.35120631e-2, 9.79340025e1, 6.96166300e-1,
        4.07942600e-1, 8.97479400e-1],        # [0] SiO2
    [4.75722038e-3, 2.37055446e-2, 1.40231330e2, 8.06866420e-1,
        7.18158480e-1, 8.54168310e-1],        # [1] GeO2
    [3.79061862e-3, 1.43810462e-2, 7.49374334e1, 6.95790000e-1,
        4.52497000e-1, 7.12513000e-1],        # [2] P2O2
    [3.83161000e-3, 1.52922902e-2, 8.27910731e1, 6.90618000e-1,
        4.01996000e-1, 8.98817000e-1],        # [3] 13.3% B2O3
    [4.65492353e-3, 1.35629316e-2, 9.98741796e1, 6.91116000e-1,
        3.99166000e-1, 8.90423000e-1],        # [4] F
    [8.90362088e-3, 8.72094500e-3, 3.59958241e1, 7.96468000e-1,
        4.97614000e-1, 3.58924000e-1],        # [5] NaO:B2O3
    [4.54079433e-2, -1.35241091e-5, 3.09549568e2, 1.30355147e-1,
        9.13764925e-1, 1.14207828],          # [6] ABCY
    [-1.55905386e-4, 7.32455962e-3, 5.96822762e2, 2.03974072e-5,
        1.25885153, 2.11857374],             # [7] HBL
    [3.31382978e-4, 9.53013988e-3, 3.85595295e2, 3.50883275e-1,
        9.36323861e-1, 1.45963548],           # [8] ZBG
    [1.49169281e-8, 8.95628044e-3, 2.39968296e2, 3.28391032e-2,
        1.25579928, 8.97176663e-1],           # [9] ZBLA
    [-2.40488039e-2, 1.73740457e-2, 4.02611805e2, 3.05900633e-1,
        9.18318740e-1, 1.50695421],          # [10] ZBLAN
    [0.004981838, 0.01375664, 97.93353, 0.6910021, 0.4022430,
        0.9439644],                             # [11] 5.2% B2O3
    [0.005202431, 0.01287730, 97.93401, 0.7058489, 0.4176021,
        0.8952753],                             # [12] 10.5% P2O2
    [6.00069867e-3, 2.00179144e-2, 103.560653, 1.03961212,
        0.231792344, 1.01046945],                  # [13] Schott N-BK7
    [4.67914826e-3, 1.35120631e-2, 97.9340025, 0.696166300,
        0.407942600, 0.897479400],                # [14] fused silica
    [5.2799261e-3, 1.42382647e-2, 325.017834, 1.43134930, 0.65054713,
        5.3414021],                     # [15] sapphire (ord. wave)
    [5.48041129e-3, 1.47994281e-2, 402.89514, 1.5039759, 0.55069141,
        6.5927379],                      # [16] sapphire (eo wave)
    [1.88217800e-3, 8.95188847e-3, 5.66135591e+2, 4.87551080e-1,
        3.98750310e-1, 2.31203530],          # [17] MgF2 (ord wave)
    [1.35737865e-3, 8.23767167e-3, 5.65107755e+2, 4.13440230e-1,
        5.04974990e-1, 2.49048620],          # [18] MgF2 (eo wave)
    [2.52642999e-3, 1.00783328e-2, 1.20055597e+3, 5.67588800e-1,
        4.71091400e-1, 3.84847230],          # [19] CaF2
    [9.97743871E-03, 4.70450767E-02, 1.11886764E+02, 1.34533359E+00,
        2.09073176E-01, 9.37357162E-01],  # Schott F2
    [9.58633048E-03, 4.57627627E-02, 1.15011883E+02, 1.31044630E+00,
        1.96034260E-01, 9.66129770E-01],  # Schott F5
    [5.20142470E-03, 1.58938446E-02, 9.59109448E+01, 9.09362180E-01,
        2.79077054E-01, 8.91813298E-01],  # Schott FK5HTi
    [8.09424251E-03, 3.86051284E-02, 1.04747730E+02, 1.15687082E+00,
        6.42625444E-02, 8.72376139E-01],  # Schott K10
    [7.20341707E-03, 2.69835916E-02, 1.00384588E+02, 1.12735550E+00,
        1.24412303E-01, 8.27100531E-01],  # Schott K7
    [1.03159999E-02, 4.69216348E-02, 8.25078509E+01, 1.66842615E+00,
        2.98512803E-01, 1.07743760E+00],  # Schott LAFN7
    [1.35670404E-02, 5.45803020E-02, 1.67904715E+02, 2.45505861E+00,
        4.53006077E-01, 2.38513080E+00],  # Schott LASF35
    [9.29854416E-03, 4.49135769E-02, 1.10493685E+02, 1.28035628E+00,
        1.63505973E-01, 8.93930112E-01],  # Schott LF5
    [9.39886260E-03, 4.52566659E-02, 1.10544829E+02, 1.28552924E+00,
        1.58357622E-01, 8.92175122E-01],  # Schott LF5HTi
    [8.57807248E-03, 4.20143003E-02, 1.07593060E+02, 1.21640125E+00,
        1.33664540E-01, 8.83399468E-01],  # Schott LLF1
    [8.70432098E-03, 4.27325235E-02, 1.08049968E+02, 1.22510445E+00,
        1.25155671E-01, 8.92236751E-01],  # Schott LLF1HTi
    [9.26681282E-03, 4.24489805E-02, 1.05613573E+02, 1.58514950E+00,
        1.43559385E-01, 1.08521269E+00],  # Schott N-BAF10
    [9.42015382E-03, 5.31087291E-02, 1.10278856E+02, 1.42056328E+00,
        1.02721269E-01, 1.14380976E+00],  # Schott N-BAF4
    [9.42734715E-03, 4.30826500E-02, 1.24889868E+02, 1.51503623E+00,
        1.53621958E-01, 1.15427909E+00],  # Schott N-BAF51
    [9.07800128E-03, 5.08212080E-02, 1.05691856E+02, 1.43903433E+00,
        9.67046052E-02, 1.09875818E+00],  # Schott N-BAF52
    [6.44742752E-03, 2.22284402E-02, 1.07297751E+02, 1.12365662E+00,
        3.09276848E-01, 8.81511957E-01],  # Schott N-BAK1
    [5.92383763E-03, 2.03828415E-02, 1.13118417E+02, 1.01662154E+00,
        3.19903051E-01, 9.37232995E-01],  # Schott N-BAK2
    [7.79980626E-03, 3.15631177E-02, 1.05965875E+02, 1.28834642E+00,
        1.32817724E-01, 9.45395373E-01],  # Schott N-BAK4
    [7.96596450E-03, 3.30672072E-02, 1.09197320E+02, 1.31004128E+00,
        1.42038259E-01, 9.64929351E-01],  # Schott N-BALF4
    [8.25815975E-03, 4.41920027E-02, 1.07097324E+02, 1.28385965E+00,
        7.19300942E-02, 1.05048927E+00],  # Schott N-BALF5
    [1.08435729E-02, 5.62278762E-02, 1.31339700E+02, 1.53652081E+00,
        1.56971102E-01, 1.30196815E+00],  # Schott N-BASF2
    [1.04485644E-02, 4.99394756E-02, 1.18961472E+02, 1.65554268E+00,
        1.71319770E-01, 1.33664448E+00],  # Schott N-BASF64
    [5.16900822E-03, 1.61190045E-02, 9.97575331E+01, 8.88308131E-01,
        3.28964475E-01, 9.84610769E-01],  # Schott N-BK10
    [9.95906143E-03, 5.46931752E-02, 1.19248346E+02, 1.39757037E+00,
        1.59201403E-01, 1.26865430E+00],  # Schott N-F2
    [4.75111955E-03, 1.49814849E-02, 9.78600293E+01, 8.44309338E-01,
        3.44147824E-01, 9.10790213E-01],  # Schott N-FK5
    [4.72301995E-03, 1.53575612E-02, 1.68681330E+02, 9.71247817E-01,
        2.16901417E-01, 9.04651666E-01],  # Schott N-FK51A
    [3.39065607E-03, 1.17551189E-02, 2.12842145E+02, 7.38042712E-01,
        3.63371967E-01, 9.89296264E-01],  # Schott N-FK58
    [6.61099503E-03, 2.41108660E-02, 1.11982777E+02, 1.08511833E+00,
        1.99562005E-01, 9.30511663E-01],  # Schott N-K5
    [8.39154696E-03, 4.04010786E-02, 1.12572446E+02, 1.19286778E+00,
        8.93346571E-02, 9.20819805E-01],  # Schott N-KF9
    [8.40298480E-03, 3.44239720E-02, 8.84310532E+01, 1.33222450E+00,
        2.89241610E-01, 1.15161734E+00],  # Schott N-KZFS11
    [7.47170505E-03, 3.08053556E-02, 7.01731084E+01, 1.23697554E+00,
        1.53569376E-01, 9.03976272E-01],  # Schott N-KZFS2
    [8.76282070E-03, 3.71767201E-02, 9.03866994E+01, 1.35055424E+00,
        1.97575506E-01, 1.09962992E+00],  # Schott N-KZFS4
    [9.86143816E-03, 4.45477583E-02, 1.06436258E+02, 1.47460789E+00,
        1.93584488E-01, 1.26589974E+00],  # Schott N-KZFS5
    [1.08808630E-02, 4.94207753E-02, 1.31009163E+02, 1.62693651E+00,
        2.43698760E-01, 1.62007141E+00],  # Schott N-KZFS8
    [1.01711622E-02, 4.42431765E-02, 1.00687748E+02, 1.80984227E+00,
        1.57295550E-01, 1.09300370E+00],  # Schott N-LAF2
    [9.33322280E-03, 3.45637762E-02, 8.32404866E+01, 1.87134529E+00,
        2.50783010E-01, 1.22048639E+00],  # Schott N-LAF21
    [9.27313493E-03, 3.58201181E-02, 8.73448712E+01, 1.79653417E+00,
        3.11577903E-01, 1.15981863E+00],  # Schott N-LAF33
    [8.72810026E-03, 2.93020832E-02, 8.51780644E+01, 1.75836958E+00,
        3.13537785E-01, 1.18925231E+00],  # Schott N-LAF34
    [7.50943203E-03, 2.60046715E-02, 8.05945159E+01, 1.51697436E+00,
        4.55875464E-01, 1.07469242E+00],  # Schott N-LAF35
    [1.07925580E-02, 5.38626639E-02, 1.06268665E+02, 1.74028764E+00,
        2.26710554E-01, 1.32525548E+00],  # Schott N-LAF7
    [8.86014635E-03, 3.63416509E-02, 8.29009069E+01, 1.72878017E+00,
        1.69257825E-01, 1.19386956E+00],  # Schott N-LAK10
    [5.77031797E-03, 2.00401678E-02, 9.54873482E+01, 1.17365704E+00,
        5.88992398E-01, 9.78014394E-01],  # Schott N-LAK12
    [7.46098727E-03, 2.42024834E-02, 8.09565165E+01, 1.50781212E+00,
        3.18866829E-01, 1.14287213E+00],  # Schott N-LAK14
    [6.02075682E-03, 1.96862889E-02, 8.84370099E+01, 1.22718116E+00,
        4.20783743E-01, 1.01284843E+00],  # Schott N-LAK21
    [5.85778594E-03, 1.98546147E-02, 1.00834017E+02, 1.14229781E+00,
        5.35138441E-01, 1.04088385E+00],  # Schott N-LAK22
    [6.70283452E-03, 2.19416210E-02, 8.07407701E+01, 1.42288601E+00,
        5.93661336E-01, 1.16135260E+00],  # Schott N-LAK33B
    [5.89278062E-03, 1.97509041E-02, 7.88894174E+01, 1.26661442E+00,
        6.65919318E-01, 1.12496120E+00],  # Schott N-LAK34
    [6.10105538E-03, 2.01388334E-02, 9.06380380E+01, 1.23679889E+00,
        4.45051837E-01, 1.01745888E+00],  # Schott N-LAK7
    [6.20023871E-03, 2.16465439E-02, 8.25827736E+01, 1.33183167E+00,
        5.46623206E-01, 1.19084015E+00],  # Schott N-LAK8
    [7.24270156E-03, 2.43353131E-02, 8.54686868E+01, 1.46231905E+00,
        3.44399589E-01, 1.15508372E+00],  # Schott N-LAK9
    [9.82060155E-03, 3.44713438E-02, 1.10739863E+02, 1.96485075E+00,
        4.75231259E-01, 1.48360109E+00],  # Schott N-LASF31A
    [1.09583310E-02, 4.74551603E-02, 9.69085286E+01, 1.98550331E+00,
        2.74057042E-01, 1.28945661E+00],  # Schott N-LASF40
    [9.10368219E-03, 3.39247268E-02, 9.33580595E+01, 1.86348331E+00,
        4.13307255E-01, 1.35784815E+00],  # Schott N-LASF41
    [1.04001413E-02, 4.47505292E-02, 8.74375690E+01, 1.93502827E+00,
        2.36629350E-01, 1.26291344E+00],  # Schott N-LASF43
    [8.72506277E-03, 3.08085023E-02, 9.27743824E+01, 1.78897105E+00,
        3.86758670E-01, 1.30506243E+00],  # Schott N-LASF44
    [1.12171920E-02, 5.05134972E-02, 1.47106505E+02, 1.87140198E+00,
        2.67777879E-01, 1.73030008E+00],  # Schott N-LASF45
    [1.23595524E-02, 5.60610282E-02, 1.07047718E+02, 2.16701566E+00,
        3.19812761E-01, 1.66004486E+00],  # Schott N-LASF46A
    [1.25805384E-02, 5.67191367E-02, 1.05316538E+02, 2.17988922E+00,
        3.06495184E-01, 1.56882437E+00],  # Schott N-LASF46B
    [1.21426017E-02, 5.38736236E-02, 1.56530829E+02, 2.00029547E+00,
        2.98926886E-01, 1.80691843E+00],  # Schott N-LASF9
    [5.85597402E-03, 1.94072416E-02, 1.40537046E+02, 1.15610775E+00,
        1.53229344E-01, 7.85618966E-01],  # Schott N-PK51
    [5.16800155E-03, 1.66658798E-02, 1.38964129E+02, 1.02960700E+00,
        1.88050600E-01, 7.36488165E-01],  # Schott N-PK52A
    [4.69824067E-03, 1.61818463E-02, 1.04374975E+02, 8.87272110E-01,
        4.89592425E-01, 1.04865296E+00],  # Schott N-PSK3
    [7.06416337E-03, 2.33251345E-02, 9.74847345E+01, 1.38121836E+00,
        1.96745645E-01, 8.86089205E-01],  # Schott N-PSK53A
    [1.19654879E-02, 5.90589722E-02, 1.35521676E+02, 1.60865158E+00,
        2.37725916E-01, 1.51530653E+00],  # Schott N-SF1
    [1.22241457E-02, 5.95736775E-02, 1.47468793E+02, 1.62153902E+00,
        2.56287842E-01, 1.64447552E+00],  # Schott N-SF10
    [1.31887070E-02, 6.23068142E-02, 1.55236290E+02, 1.73759695E+00,
        3.13747346E-01, 1.89878101E+00],  # Schott N-SF11
    [1.30512113E-02, 6.13691880E-02, 1.49517689E+02, 1.69022361E+00,
        2.88870052E-01, 1.70451870E+00],  # Schott N-SF14
    [1.16507014E-02, 5.97856897E-02, 1.32709339E+02, 1.57055634E+00,
        2.18987094E-01, 1.50824017E+00],  # Schott N-SF15
    [1.09019098E-02, 5.85683687E-02, 1.27404933E+02, 1.47343127E+00,
        1.63681849E-01, 1.36920899E+00],  # Schott N-SF2
    [1.26793450E-02, 6.02038419E-02, 1.45760496E+02, 1.67780282E+00,
        2.82849893E-01, 1.63539276E+00],  # Schott N-SF4
    [1.12547560E-02, 5.88995392E-02, 1.29141675E+02, 1.52481889E+00,
        1.87085527E-01, 1.42729015E+00],  # Schott N-SF5
    [1.41749518E-02, 6.40509927E-02, 1.77389795E+02, 1.87543831E+00,
        3.73757490E-01, 2.30001797E+00],  # Schott N-SF57
    [1.33714182E-02, 6.17533621E-02, 1.74017590E+02, 1.77931763E+00,
        3.38149866E-01, 2.08734474E+00],  # Schott N-SF6
    [1.47053225E-02, 6.92998276E-02, 1.61817601E+02, 2.02459760E+00,
        4.70187196E-01, 2.59970433E+00],  # Schott N-SF66
    [1.33714182E-02, 6.17533621E-02, 1.74017590E+02, 1.77931763E+00,
        3.38149866E-01, 2.08734474E+00],  # Schott N-SF6HT
    [1.14338344E-02, 5.82725652E-02, 1.33241650E+02, 1.55075812E+00,
        2.09816918E-01, 1.46205491E+00],  # Schott N-SF8
    [6.80282081E-03, 2.19737205E-02, 1.01513232E+02, 1.17963631E+00,
        2.29817295E-01, 9.35789652E-01],  # Schott N-SK11
    [4.61716525E-03, 1.68859270E-02, 1.03736265E+02, 9.36155374E-01,
        5.94052018E-01, 1.04374583E+00],  # Schott N-SK14
    [7.04687339E-03, 2.29005000E-02, 9.27508526E+01, 1.34317774E+00,
        2.41144399E-01, 9.94317969E-01],  # Schott N-SK16
    [7.27191640E-03, 2.42823527E-02, 1.10377773E+02, 1.28189012E+00,
        2.57738258E-01, 9.68186040E-01],  # Schott N-SK2
    [7.16874107E-03, 2.46455892E-02, 1.00886364E+02, 1.32993741E+00,
        2.28542996E-01, 9.88465211E-01],  # Schott N-SK4
    [5.22730467E-03, 1.72733646E-02, 9.83594579E+01, 9.91463823E-01,
        4.95982121E-01, 9.87393925E-01],  # Schott N-SK5
    [8.23982975E-03, 3.33736841E-02, 1.06870822E+02, 1.43060270E+00,
        1.53150554E-01, 1.01390904E+00],  # Schott N-SSK2
    [9.20284626E-03, 4.23530072E-02, 1.06927374E+02, 1.59222659E+00,
        1.03520774E-01, 1.05174016E+00],  # Schott N-SSK5
    [8.69310149E-03, 4.21566593E-02, 1.11300666E+02, 1.44857867E+00,
        1.17965926E-01, 1.06937528E+00],  # Schott N-SSK8
    [6.76601657E-03, 2.30642817E-02, 8.90498778E+01, 1.07715032E+00,
        1.68079109E-01, 8.51889892E-01],  # Schott N-ZK7
    [6.76601657E-03, 2.30642817E-02, 8.90498778E+01, 1.07509891E+00,
        1.68895044E-01, 8.60503983E-01],  # Schott N-ZK7A
    [7.22141956E-03, 2.68216805E-02, 1.01702362E+02, 1.18318503E+00,
        8.71756426E-02, 1.03133701E+00],  # Schott P-BK7
    [9.38006396E-03, 3.60537464E-02, 8.64324693E+01, 1.76003244E+00,
        2.48286745E-01, 1.15935122E+00],  # Schott P-LAF37
    [7.15959695E-03, 2.33637446E-02, 8.83284426E+01, 1.39324260E+00,
        4.18882766E-01, 1.04380700E+00],  # Schott P-LAK35
    [1.00328203E-02, 3.87095168E-02, 9.45421507E+01, 1.85543101E+00,
        3.15854649E-01, 1.28561839E+00],  # Schott P-LASF47
    [9.99234757E-03, 3.87437988E-02, 9.58967681E+01, 1.84910553E+00,
        3.29828674E-01, 1.30400901E+00],  # Schott P-LASF50
    [9.88495571E-03, 3.78097402E-02, 9.78415430E+01, 1.84568806E+00,
        3.39001600E-01, 1.32418921E+00],  # Schott P-LASF51
    [1.68838419E-02, 7.16086325E-02, 1.18707479E+02, 2.33300670E+00,
        4.52961396E-01, 1.25172339E+00],  # Schott P-SF68
    [1.21696677E-02, 6.00710405E-02, 1.45651908E+02, 1.62594647E+00,
        2.35927609E-01, 1.67434623E+00],  # Schott P-SF69
    [1.16582670E-02, 5.82087757E-02, 1.30748028E+02, 1.55370411E+00,
        2.06332561E-01, 1.39708831E+00],  # Schott P-SF8
    [7.40877235E-03, 2.54563489E-02, 1.07751087E+02, 1.31053414E+00,
        1.69376189E-01, 1.10987714E+00],  # Schott P-SK57
    [7.36408831E-03, 2.55786047E-02, 1.06726060E+02, 1.30536483E+00,
        1.71434328E-01, 1.10117219E+00],  # Schott P-SK57Q1
    [7.20717498E-03, 2.45659595E-02, 1.02739728E+02, 1.31678410E+00,
        1.71154756E-01, 1.12501473E+00],  # Schott P-SK58A
    [7.84382378E-03, 2.87769365E-02, 1.05373397E+02, 1.40790442E+00,
        1.43381417E-01, 1.16513947E+00],  # Schott P-SK60
    [1.21481001E-02, 5.34549042E-02, 1.12174809E+02, 1.55912923E+00,
        2.84246288E-01, 9.68842926E-01],  # Schott SF1
    [1.27534559E-02, 5.81983954E-02, 1.16607680E+02, 1.61625977E+00,
        2.59229334E-01, 1.07762317E+00],  # Schott SF10
    [1.36068604E-02, 6.15960463E-02, 1.21922711E+02, 1.73848403E+00,
        3.11168974E-01, 1.17490871E+00],  # Schott SF11
    [1.05795466E-02, 4.93226978E-02, 1.12405955E+02, 1.40301821E+00,
        2.31767504E-01, 9.39056586E-01],  # Schott SF2
    [1.25502104E-02, 5.44559822E-02, 1.17652222E+02, 1.61957826E+00,
        3.39493189E-01, 1.02566931E+00],  # Schott SF4
    [1.11826126E-02, 5.08594669E-02, 1.12041888E+02, 1.46141885E+00,
        2.47713019E-01, 9.49995832E-01],  # Schott SF5
    [1.33874699E-02, 5.79561608E-02, 1.21616024E+02, 1.70579259E+00,
        3.44223052E-01, 1.09601828E+00],  # Schott SF56A
    [1.43704198E-02, 5.92801172E-02, 1.21419942E+02, 1.81651371E+00,
        4.28893641E-01, 1.07186278E+00]  # Schott SF57
]

ALL_GLASS_NAMES = np.array([
    "SiO2", "GeO2", "9.1% P2O2", "13.3% B2O3", "1.0% F",
    "16.9% Na2O : 32.5% B2O3", "ABCY", "HBL", "ZBG", "ZBLA", "ZBLAN",
    "5.2% B2O3", "10.5% P2O2", "N-BK7", "fused silica",
    "sapphire (ordinary)", "sapphire (extraordinary)", "MgF2 (ordinary)",
    "MgF2 (extraordinary)", "CaF2", "F2", "F5", "FK5HTi", "K10", "K7", "LAFN7",
    "LASF35", "LF5", "LF5HTi", "LLF1", "LLF1HT", "N-BAF1", "N-BAF4", "N-BAF5",
    "N-BAF5", "N-BAK1", "N-BAK2", "N-BAK4", "N-BALF", "N-BALF", "N-BASF", "N-BASF",
    "N-BK10", "N-F2", "N-FK5", "N-FK51", "N-FK58", "N-K5", "N-KF9", "N-KZFS",
    "N-KZFS", "N-KZFS", "N-KZFS", "N-KZFS", "N-LAF2", "N-LAF2", "N-LAF3",
    "N-LAF3", "N-LAF3", "N-LAF7", "N-LAK1", "N-LAK1", "N-LAK1", "N-LAK2",
    "N-LAK2", "N-LAK3", "N-LAK3", "N-LAK7", "N-LAK8", "N-LAK9", "N-LASF",
    "N-LASF", "N-LASF", "N-LASF", "N-LASF", "N-LASF", "N-LASF", "N-LASF",
    "N-LASF", "N-PK51", "N-PK52", "N-PSK3", "N-PSK5", "N-SF1", "N-SF10",
    "N-SF11", "N-SF14", "N-SF15", "N-SF2", "N-SF4", "N-SF5", "N-SF57",
    "N-SF6", "N-SF66", "N-SF6H", "N-SF8", "N-SK11", "N-SK14", "N-SK16",
    "N-SK2", "N-SK4", "N-SK5", "N-SSK2", "N-SSK5", "N-SSK8", "N-ZK7",
    "N-ZK7A", "P-BK7", "P-LAF3", "P-LAK3", "P-LASF", "P-LASF", "P-LASF",
    "P-SF68", "P-SF69", "P-SF8", "P-SK57", "P-SK57", "P-SK58", "P-SK60",
    "SF1", "SF10", "SF11", "SF2", "SF4", "SF5", "SF56A", "SF57"])


def glass(i):
    """
    returns an array of Sellmeier coefficients for glass with index i
    Use like this
        lambdas = np.linspace(1000,1700,50)*1e-9 # [m]
        i = ofiber.glass_index("SiO2")
        glass = ofiber.refraction.glass(i)       # SiO2
        n = ofiber.refraction.n(glass,lambdas)
        plt.plot(lambdas*1e9, n)
    """
    return _GLASS[i]


def glass_name(i):
    """
    Look up the name of the glass with index i

    (A list of all possible names is in the array
        ofiber.refraction.ALL_GLASS_NAMES
    )
    """
    return ALL_GLASS_NAMES[i]


def find_glass(name):
    """
    Look up the index of the glass with a particular name

    The index of the first glass that matches the string is returned.
    Matching is case insensitive

    Arg:
        name: string containing
    Returns:
        the index of the first glass that matches
    """
    target = name.upper()

    for i, i_glass_name in enumerate(ALL_GLASS_NAMES):
        if target in i_glass_name.upper():
            return i

    print("'%s' not found in " % target, ALL_GLASS_NAMES)
    return 0


def doped_glass(x):
    """
    Calculate Sellmeier coefficients for SiO_2 doped with GeO_2

    Arg:
        x = molar fraction of GeO_2 in the system,
                x * GeO_2 : (1 - x) * SiO_2
    Returns:
        Sellmeier coefficients for doped glass (array of six values)
    """
    SA = np.array([0.6961663, 0.4079426, 0.8974794])
    SL = np.array([0.0684043, 0.1162414, 9.896161])
    GA = np.array([0.80686642, 0.71815848, 0.85416831])
    GL = np.array([0.068972606, 0.15396605, 11.841931])
    a = (SL + x * (GL - SL))**2
    b = abs(SA + x * (GA - SA))
    return np.concatenate([a, b])


def doped_glass_name(x):
    """
    Create a string the name describing the GeO2 doped glass

    Arg:
        x = molar fraction of GeO_2 in the system,
    Returns:
        string describing the doped glass
    """
    if x == 0:
        return r'SiO$_2$'
    if x == 1:
        return r'GeO$_2$'

    return r'%.2f GeO$_2$ : %.2f SiO$_2$' % (x, 1 - x)


def _sellmeier(b, c, lambda0):
    """
    Calculates the index of refraction using the Sellmeier equation

    This is intended as a private method.

    Args:
        b : array of three Sellmeier Coefficients  [--]
        c : array of three Sellmeier Coefficients  [microns**2]
        lambda0 : wavelength in vacuum             [m]
    Returns:
        returns the index of refraction at lambda0 [-]
    """
    lam2 = lambda0**2 * 1e12  # um**2
    nsq = 1
    for i in range(3):
        nsq += b[i] * lam2 / (lam2 - c[i])

    return np.sqrt(nsq)


def _d_sellmeier(b, c, lambda0):
    """
    Calculates the first derivative (wrt wavelength) of the Sellmeier equation

    This is a private method.

    Args:
        b : array of three Sellmeier Coefficients  [--]
        c : array of three Sellmeier Coefficients  [microns**2]
        lambda0 : wavelength in vacuum             [m]
    Returns:
        returns the first derivative of the index of refraction at lambda0 [1/m]
    """
    n1 = _sellmeier(b, c, lambda0)
    lam = lambda0 * 1e6  # microns
    lam2 = lam**2

    dy = 0
    for i in range(3):
        dy -= b[i] * c[i] / (lam2 - c[i])**2  # 1/um**2

    dy *= lam / n1                            # 1/um
    dy *= 1e6                                 # 1/m

    return dy


def _d2_sellmeier(b, c, lambda0):
    """
    Calculates the second derivative (wrt wavelength) of the Sellmeier equation

    This is a private method.

    Args:
        b : array of three Sellmeier Coefficients  [--]
        c : array of three Sellmeier Coefficients  [microns**2]
        lambda0 : wavelength in vacuum             [m]
    Returns:
        returns the second derivative of the refractive index at lambda0 [1/m]
    """
    nn = _sellmeier(b, c, lambda0)    # index of refraction
    lam = lambda0 * 1e6               # needed because Sellmeier uses [um]
    lam2 = lam**2                     # [um**2]

    dy = 0
    d2y = 0
    for i in range(3):
        dy = b[i] * c[i] / (lam2 - c[i])**2                        # 1/um
        d2y += b[i] * c[i] * (3 * lam2 + c[i]) / (lam2 - c[i])**3  # 1/um**2

    total = d2y / nn - lam2 * dy**2 / nn**3                        # 1/um**2
    total *= 1e12                                                  # 1/m**2

    return total


def n(glass_coef, lambda0):
    """
    Calculates index of refraction for a glass at lambda0

    Args:
        glass_coef: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]
    Returns:
        index of refraction [--]
    """
    return _sellmeier(glass_coef[3:6], glass_coef[0:3], lambda0)


def dn(glass_coef, lambda0):
    """
    Calculates the first derivative of the refractiv index w.r.t. wavelength

    Args:
        glass_coef: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]
    Returns:
        the first derivative of index of refraction [1/m]
    """
    return _d_sellmeier(glass_coef[3:6], glass_coef[0:3], lambda0)


def d2n(glass_coef, lambda0):
    """
    Calculates the second derivative of the refractiv index w.r.t. wavelength

    Args:
        glass_coef: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]
    Returns:
        the second derivative of index of refraction [1/m**2]
    """
    return _d2_sellmeier(glass_coef[3:6], glass_coef[0:3], lambda0)


def n_group(glass_coef, lambda0):
    """
    Calculates group index of refraction at lambda0

    Args:
        glass: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]
    Returns:
        group index of refraction [--]
    """
    return n(glass_coef, lambda0) - lambda0*dn(glass_coef, lambda0)


def n_air(lambda0, temperature=15):
    """
    Calculates refractive index of air at atmospheric pressure, following
    equation given on page 4 of Smith, Modern Optical Engineering

    Args:
        lambda0: wavelength in vacuum [m]
        temperature: degrees celsius
    Returns:
        index of refraction [--]
    """
    nu = 1/(lambda0 * 1e6)
    n15 = 1e-8 * (8342.1 + 2406030/(130 - nu**2) + 15996/(38.9 - nu**2))
    if temperature == 15:
        return 1 + n15

    return 1 + 1.0549*n15/(1+0.00366*temperature)

# code used to generate the Sellmeier coefficients for other_glass below
# data straight from fleming 1978
# flemingzz=np.array([
#     [0.696750, 0.069066, 0.408218, 0.115662, 0.890815, 9.900559],
#     [0.711040, 0.064270, 0.451885, 0.129408, 0.704048, 9.425478],
#     [0.695790, 0.061568, 0.452497, 0.119921, 0.712513, 8.656641],
#     [0.690618, 0.061900, 0.401996, 0.123662, 0.898817, 9.098960],
#     [0.691116, 0.068227, 0.399166, 0.116460, 0.890423, 9.993707],
#     [0.796468, 0.094359, 0.497614, 0.093386, 0.358924, 5.999652]]).T
#
# rearrange so each row is a type of glass with Sellmeier coefficients needed
# b1,l1,b2,l2,b3,l3=flemingzz
# a1=l1**2
# a2=l2**2
# a3=l3**2
# fleming_all=np.array([a1,a2,a3,b1,b2,b3]).T
# fleming_names = np.array(["Quenched SiO$_2$","13.5% Ge$O_2$",\
#     "9.1% P$_2$O$_2$","13.3% B$_2$O$_3$","1.0% F",\
#     "16.9% Na$_2$O : 32.5% B$_2$O$_3$"])
#
# extract glasses that are not SiO2:GeO2 mixtures
# glass = fleming_all[2:6]
# glass_names = fleming_names[2:6]


# data straight from Ghatak
# pure_sio2 = [0.004679148,0.01351206,97.93400,0.6961663,0.4079426,0.8974794]
# geo2_63 = [0.007290464,0.01050294,97.93428,0.7083952,0.4203993,0.8663412]
# geo2_193 = [0.005847345,0.01552717,97.93484,0.7347008,0.4461191,0.8081698]
# b203 = [0.004981838,0.01375664,97.93353,0.6910021,0.4022430,0.9439644]
# p2o3 = [0.005202431,0.01287730,97.93401,0.7058489,0.4176021,0.8952753]


# def paek_refraction(lambda0):
#     """
#     Return the index of refraction of silica fiber at
#     vacuum wavelength lambda0 [m]
#     result is dimensionless
#     """
#     ell = 0.035       # um**2
#     c0 = 1.4508554
#     c1 = -0.0031268
#     c2 = -0.0000381
#     c3 = 0.0030270
#     c4 = -0.0000779
#     c5 = 0.0000018
#     lam = lambda0 * 1e6   # now in microns
#     den = lam**2 - ell
#     n = c0 + c1 * lam**2 + c2 * lam**4 + c3 / den + c4 / den**2 + c5 / den**3
#     return n
#
#
# def d_paek_refraction(lambda0):
#     """
#     Return the first derivative of index of refraction of silica fiber at
#     vacuum wavelength lambda0 [m] with respect to wavelength
#     result has units of [1/m]
#     """
#     ell = 0.035
#     c1 = -0.0031268
#     c2 = -0.0000381
#     c3 = 0.0030270
#     c4 = -0.0000779
#     c5 = 0.0000018
#     lam = lambda0 * 1e6   # microns
#     den = lam**2 - ell
#     n = 2 * c1 * lam + 4 * c2 * lam**3 - 2 * c3 * lam / den**2
#     n += - 4 * c4 * lam / den**3 - 6 * c5 * lam / den**4      # um**-1
#     n *= 1e6                                                  # m**-1
#     return n
#
#
# def d2_paek_refraction(lambda0):
#     """
#     Return the second derivative of index of refraction of silica fiber at
#     vacuum wavelength lambda0 [m] with respect to wavelength
#     result has units of [1/m**2]
#     """
#     ell = 0.035
#     c1 = -0.0031268
#     c2 = -0.0000381
#     c3 = 0.0030270
#     c4 = -0.0000779
#     c5 = 0.0000018
#     lam = lambda0 * 1e6   # microns
#     den = lam**2 - ell
#     n = 2 * c1 + 12 * c2 * lam**2 - 2 * c3 / \
#         den**2 - 4 * c4 / den**3 - 6 * c5 / den**4
#     n += 48 * c5 * lam**2 / den**5 + 24 * c4 * \
#         lam**2 / den**4 + 8 * c3 * lam**2 / den**3
#     n *= 1e12           # m**-2
#     return n
#
# ### method for converting from mendez coefficients to sellmeier coefficients
# ### just a simple curve fit
#
# from scipy.optimize import curve_fit
#
# def sell(x,a1,a2,a3,b1,b2,b3):
#     x2 = x**2
#     return np.sqrt(1 + b1*x2/(x2 - a1) + b2*x2/(x2 - a2) + b3*x2/(x2 - a3))
#
#
# lambda0 = np.linspace(1500,3500,50)*1e-9 # [m]
# for i in range(len(ofr.mendez_glass)):
#     glass = ofr.mendez_glass[i]
#     n = ofr.mendez_refraction(glass,lambda0)
#     popt, pcov = curve_fit(sell, lambda0*1e6, n, p0=ofr.sell_glass[0])
#     plt.plot(lambda0*1e9, n-sell(lambda0*1e6, *popt), 'r-', label='fit')
#     print(popt)

# mendez_glass =[
# [7.67742e-6, 2.16195e-3, 1.42969, -1.28304e-3, -5.35487e-6],
# [-28.61020e-6, 3.11470e-3, 1.50294, -1.17821e-3, -2.64123e-6],
# [93.67070e-6, 2.94329e-3, 1.51236, -1.25045e-3, -4.01026e-6],
# [-300.80370e-6, 4.03214e-3, 1.51272, -1.21921e-3, -6.77630e-6],
# [93.67070e-6, 2.94329e-3, 1.49136, -1.25045e-3, -4.01026e-6]
# ]
# mendezall_glass_namess=["ABCY", "HBL", "ZBG", "ZBLA", "ZBLAN"]
#
# def mendez_refraction(glass, lambda0):
#     """
#     returns the index of refraction using the Mendez equation
#     lambda0 is in [m]
#     """
#     lam = lambda0 * 1e6   # microns
#     sum = glass[0] * lam**-4
#     sum += glass[1] * lam**-2
#     sum += glass[2]
#     sum += glass[3] * lam**2
#     sum += glass[4] * lam**4
#     return sum
#
#
# def d_mendez_refraction(glass, lambda0):
#     """
#     returns first derivative of the  index of refraction using the Mendez eqn
#     lambda0 is in [m]
#     """
#     sum = 0
#     lam = lambda0 * 1e6   # microns
#     sum = -4 * glass[0] * lam**-5
#     sum += -2 * glass[1] * lam**-3
#     sum += 2 * glass[3] * lam
#     sum += 4 * glass[4] * lam**3
#     sum *= 1e6  # [1/m]
#     return sum
#
#
# def d2_mendez_refraction(glass, lambda0):
#     """
#     returns second derivative of the  index of refraction using the Mendez eqn
#     lambda0 is in [m]
#     """
#     sum = 0
#     lam = lambda0 * 1e6   # microns
#     sum = 20 * glass[0] * lam**-6
#     sum += 6 * glass[1] * lam**-4
#     sum += 2 * glass[3]
#     sum += 12 * glass[4] * lam**2
#     sum *= 1e12  # [1/m]
#     return sum
