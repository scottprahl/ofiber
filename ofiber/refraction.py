# pylint: disable=invalid-name
# pylint: disable=line-too-long
# pylint: disable=consider-using-f-string
"""
Useful routines for index of refraction based on Sellmeier coefficients.

See <https://ofiber.readthedocs.io> for usage examples.

Included are Sellmeier coefficients for a bunch of common glasses used
in optical fibers.  Functions are provided to generate the index of refraction,
as well as its first and second derivatives.

Use like this::
    import ofiber as of
    lambdas = np.linspace(1000,1700,50)*1e-9 # [m]
    i = of.find_glass("SiO2")
    glass = of.glass(i)
    n = of.n(glass,lambdas)
    plt.plot(lambdas*1e9, n)
"""

import numpy as np


__all__ = (
    "ALL_GLASS_NAMES",
    "d2n",
    "dn",
    "doped_glass",
    "doped_glass_name",
    "glass",
    "glass_name",
    "find_glass",
    "n",
    "n_group",
    "n_air",
)

_GLASS = [
    # format is [c1, c2, c3, b1, b2, b3]
    # where b1, b2, b3 are unitless and c1, c2, c3 have units of [microns**2]
    [4.67914826e-3, 1.35120631e-2, 9.79340025e1, 6.96166300e-1, 4.07942600e-1, 8.97479400e-1],  # [0] SiO₂
    [4.75722038e-3, 2.37055446e-2, 1.40231330e2, 8.06866420e-1, 7.18158480e-1, 8.54168310e-1],  # [1] GeO₂
    [3.79061862e-3, 1.43810462e-2, 7.49374334e1, 6.95790000e-1, 4.52497000e-1, 7.12513000e-1],  # [2] P2O₂
    [3.83161000e-3, 1.52922902e-2, 8.27910731e1, 6.90618000e-1, 4.01996000e-1, 8.98817000e-1],  # [3] 13.3% B2O3
    [4.65492353e-3, 1.35629316e-2, 9.98741796e1, 6.91116000e-1, 3.99166000e-1, 8.90423000e-1],  # [4] F
    [8.90362088e-3, 8.72094500e-3, 3.59958241e1, 7.96468000e-1, 4.97614000e-1, 3.58924000e-1],  # [5] NaO:B2O3
    [4.54079433e-2, -1.35241091e-5, 3.09549568e2, 1.30355147e-1, 9.13764925e-1, 1.14207828],  # [6] ABCY
    [-1.55905386e-4, 7.32455962e-3, 5.96822762e2, 2.03974072e-5, 1.25885153, 2.11857374],  # [7] HBL
    [3.31382978e-4, 9.53013988e-3, 3.85595295e2, 3.50883275e-1, 9.36323861e-1, 1.45963548],  # [8] ZBG
    [1.49169281e-8, 8.95628044e-3, 2.39968296e2, 3.28391032e-2, 1.25579928, 8.97176663e-1],  # [9] ZBLA
    [-2.40488039e-2, 1.73740457e-2, 4.02611805e2, 3.05900633e-1, 9.18318740e-1, 1.50695421],  # [10] ZBLAN
    [0.004981838, 0.01375664, 97.93353, 0.6910021, 0.4022430, 0.9439644],  # [11] 5.2% B2O3
    [0.005202431, 0.01287730, 97.93401, 0.7058489, 0.4176021, 0.8952753],  # [12] 10.5% P2O₂
    [6.00069867e-3, 2.00179144e-2, 103.560653, 1.03961212, 0.231792344, 1.01046945],  # [13] Schott N-BK7
    [4.67914826e-3, 1.35120631e-2, 97.9340025, 0.696166300, 0.407942600, 0.897479400],  # [14] fused silica
    [5.2799261e-3, 1.42382647e-2, 325.017834, 1.43134930, 0.65054713, 5.3414021],  # [15] sapphire (ord. wave)
    [5.48041129e-3, 1.47994281e-2, 402.89514, 1.5039759, 0.55069141, 6.5927379],  # [16] sapphire (eo wave)
    [1.88217800e-3, 8.95188847e-3, 5.66135591e2, 4.87551080e-1, 3.98750310e-1, 2.31203530],  # [17] MgF2 (ord wave)
    [1.35737865e-3, 8.23767167e-3, 5.65107755e2, 4.13440230e-1, 5.04974990e-1, 2.49048620],  # [18] MgF2 (eo wave)
    [2.52642999e-3, 1.00783328e-2, 1.20055597e3, 5.67588800e-1, 4.71091400e-1, 3.84847230],  # [19] CaF2
    [9.97743871e-03, 4.70450767e-02, 1.11886764e02, 1.34533359e00, 2.09073176e-01, 9.37357162e-01],  # Schott F2
    [9.58633048e-03, 4.57627627e-02, 1.15011883e02, 1.31044630e00, 1.96034260e-01, 9.66129770e-01],  # Schott F5
    [5.20142470e-03, 1.58938446e-02, 9.59109448e01, 9.09362180e-01, 2.79077054e-01, 8.91813298e-01],  # Schott FK5HTi
    [8.09424251e-03, 3.86051284e-02, 1.04747730e02, 1.15687082e00, 6.42625444e-02, 8.72376139e-01],  # Schott K10
    [7.20341707e-03, 2.69835916e-02, 1.00384588e02, 1.12735550e00, 1.24412303e-01, 8.27100531e-01],  # Schott K7
    [1.03159999e-02, 4.69216348e-02, 8.25078509e01, 1.66842615e00, 2.98512803e-01, 1.07743760e00],  # Schott LAFN7
    [1.35670404e-02, 5.45803020e-02, 1.67904715e02, 2.45505861e00, 4.53006077e-01, 2.38513080e00],  # Schott LASF35
    [9.29854416e-03, 4.49135769e-02, 1.10493685e02, 1.28035628e00, 1.63505973e-01, 8.93930112e-01],  # Schott LF5
    [9.39886260e-03, 4.52566659e-02, 1.10544829e02, 1.28552924e00, 1.58357622e-01, 8.92175122e-01],  # Schott LF5HTi
    [8.57807248e-03, 4.20143003e-02, 1.07593060e02, 1.21640125e00, 1.33664540e-01, 8.83399468e-01],  # Schott LLF1
    [8.70432098e-03, 4.27325235e-02, 1.08049968e02, 1.22510445e00, 1.25155671e-01, 8.92236751e-01],  # Schott LLF1HTi
    [9.26681282e-03, 4.24489805e-02, 1.05613573e02, 1.58514950e00, 1.43559385e-01, 1.08521269e00],  # Schott N-BAF10
    [9.42015382e-03, 5.31087291e-02, 1.10278856e02, 1.42056328e00, 1.02721269e-01, 1.14380976e00],  # Schott N-BAF4
    [9.42734715e-03, 4.30826500e-02, 1.24889868e02, 1.51503623e00, 1.53621958e-01, 1.15427909e00],  # Schott N-BAF51
    [9.07800128e-03, 5.08212080e-02, 1.05691856e02, 1.43903433e00, 9.67046052e-02, 1.09875818e00],  # Schott N-BAF52
    [6.44742752e-03, 2.22284402e-02, 1.07297751e02, 1.12365662e00, 3.09276848e-01, 8.81511957e-01],  # Schott N-BAK1
    [5.92383763e-03, 2.03828415e-02, 1.13118417e02, 1.01662154e00, 3.19903051e-01, 9.37232995e-01],  # Schott N-BAK2
    [7.79980626e-03, 3.15631177e-02, 1.05965875e02, 1.28834642e00, 1.32817724e-01, 9.45395373e-01],  # Schott N-BAK4
    [7.96596450e-03, 3.30672072e-02, 1.09197320e02, 1.31004128e00, 1.42038259e-01, 9.64929351e-01],  # Schott N-BALF4
    [8.25815975e-03, 4.41920027e-02, 1.07097324e02, 1.28385965e00, 7.19300942e-02, 1.05048927e00],  # Schott N-BALF5
    [1.08435729e-02, 5.62278762e-02, 1.31339700e02, 1.53652081e00, 1.56971102e-01, 1.30196815e00],  # Schott N-BASF2
    [1.04485644e-02, 4.99394756e-02, 1.18961472e02, 1.65554268e00, 1.71319770e-01, 1.33664448e00],  # Schott N-BASF64
    [5.16900822e-03, 1.61190045e-02, 9.97575331e01, 8.88308131e-01, 3.28964475e-01, 9.84610769e-01],  # Schott N-BK10
    [9.95906143e-03, 5.46931752e-02, 1.19248346e02, 1.39757037e00, 1.59201403e-01, 1.26865430e00],  # Schott N-F2
    [4.75111955e-03, 1.49814849e-02, 9.78600293e01, 8.44309338e-01, 3.44147824e-01, 9.10790213e-01],  # Schott N-FK5
    [4.72301995e-03, 1.53575612e-02, 1.68681330e02, 9.71247817e-01, 2.16901417e-01, 9.04651666e-01],  # Schott N-FK51A
    [3.39065607e-03, 1.17551189e-02, 2.12842145e02, 7.38042712e-01, 3.63371967e-01, 9.89296264e-01],  # Schott N-FK58
    [6.61099503e-03, 2.41108660e-02, 1.11982777e02, 1.08511833e00, 1.99562005e-01, 9.30511663e-01],  # Schott N-K5
    [8.39154696e-03, 4.04010786e-02, 1.12572446e02, 1.19286778e00, 8.93346571e-02, 9.20819805e-01],  # Schott N-KF9
    [8.40298480e-03, 3.44239720e-02, 8.84310532e01, 1.33222450e00, 2.89241610e-01, 1.15161734e00],  # Schott N-KZFS11
    [7.47170505e-03, 3.08053556e-02, 7.01731084e01, 1.23697554e00, 1.53569376e-01, 9.03976272e-01],  # Schott N-KZFS2
    [8.76282070e-03, 3.71767201e-02, 9.03866994e01, 1.35055424e00, 1.97575506e-01, 1.09962992e00],  # Schott N-KZFS4
    [9.86143816e-03, 4.45477583e-02, 1.06436258e02, 1.47460789e00, 1.93584488e-01, 1.26589974e00],  # Schott N-KZFS5
    [1.08808630e-02, 4.94207753e-02, 1.31009163e02, 1.62693651e00, 2.43698760e-01, 1.62007141e00],  # Schott N-KZFS8
    [1.01711622e-02, 4.42431765e-02, 1.00687748e02, 1.80984227e00, 1.57295550e-01, 1.09300370e00],  # Schott N-LAF2
    [9.33322280e-03, 3.45637762e-02, 8.32404866e01, 1.87134529e00, 2.50783010e-01, 1.22048639e00],  # Schott N-LAF21
    [9.27313493e-03, 3.58201181e-02, 8.73448712e01, 1.79653417e00, 3.11577903e-01, 1.15981863e00],  # Schott N-LAF33
    [8.72810026e-03, 2.93020832e-02, 8.51780644e01, 1.75836958e00, 3.13537785e-01, 1.18925231e00],  # Schott N-LAF34
    [7.50943203e-03, 2.60046715e-02, 8.05945159e01, 1.51697436e00, 4.55875464e-01, 1.07469242e00],  # Schott N-LAF35
    [1.07925580e-02, 5.38626639e-02, 1.06268665e02, 1.74028764e00, 2.26710554e-01, 1.32525548e00],  # Schott N-LAF7
    [8.86014635e-03, 3.63416509e-02, 8.29009069e01, 1.72878017e00, 1.69257825e-01, 1.19386956e00],  # Schott N-LAK10
    [5.77031797e-03, 2.00401678e-02, 9.54873482e01, 1.17365704e00, 5.88992398e-01, 9.78014394e-01],  # Schott N-LAK12
    [7.46098727e-03, 2.42024834e-02, 8.09565165e01, 1.50781212e00, 3.18866829e-01, 1.14287213e00],  # Schott N-LAK14
    [6.02075682e-03, 1.96862889e-02, 8.84370099e01, 1.22718116e00, 4.20783743e-01, 1.01284843e00],  # Schott N-LAK21
    [5.85778594e-03, 1.98546147e-02, 1.00834017e02, 1.14229781e00, 5.35138441e-01, 1.04088385e00],  # Schott N-LAK22
    [6.70283452e-03, 2.19416210e-02, 8.07407701e01, 1.42288601e00, 5.93661336e-01, 1.16135260e00],  # Schott N-LAK33B
    [5.89278062e-03, 1.97509041e-02, 7.88894174e01, 1.26661442e00, 6.65919318e-01, 1.12496120e00],  # Schott N-LAK34
    [6.10105538e-03, 2.01388334e-02, 9.06380380e01, 1.23679889e00, 4.45051837e-01, 1.01745888e00],  # Schott N-LAK7
    [6.20023871e-03, 2.16465439e-02, 8.25827736e01, 1.33183167e00, 5.46623206e-01, 1.19084015e00],  # Schott N-LAK8
    [7.24270156e-03, 2.43353131e-02, 8.54686868e01, 1.46231905e00, 3.44399589e-01, 1.15508372e00],  # Schott N-LAK9
    [9.82060155e-03, 3.44713438e-02, 1.10739863e02, 1.96485075e00, 4.75231259e-01, 1.48360109e00],  # Schott N-LASF31A
    [1.09583310e-02, 4.74551603e-02, 9.69085286e01, 1.98550331e00, 2.74057042e-01, 1.28945661e00],  # Schott N-LASF40
    [9.10368219e-03, 3.39247268e-02, 9.33580595e01, 1.86348331e00, 4.13307255e-01, 1.35784815e00],  # Schott N-LASF41
    [1.04001413e-02, 4.47505292e-02, 8.74375690e01, 1.93502827e00, 2.36629350e-01, 1.26291344e00],  # Schott N-LASF43
    [8.72506277e-03, 3.08085023e-02, 9.27743824e01, 1.78897105e00, 3.86758670e-01, 1.30506243e00],  # Schott N-LASF44
    [1.12171920e-02, 5.05134972e-02, 1.47106505e02, 1.87140198e00, 2.67777879e-01, 1.73030008e00],  # Schott N-LASF45
    [1.23595524e-02, 5.60610282e-02, 1.07047718e02, 2.16701566e00, 3.19812761e-01, 1.66004486e00],  # Schott N-LASF46A
    [1.25805384e-02, 5.67191367e-02, 1.05316538e02, 2.17988922e00, 3.06495184e-01, 1.56882437e00],  # Schott N-LASF46B
    [1.21426017e-02, 5.38736236e-02, 1.56530829e02, 2.00029547e00, 2.98926886e-01, 1.80691843e00],  # Schott N-LASF9
    [5.85597402e-03, 1.94072416e-02, 1.40537046e02, 1.15610775e00, 1.53229344e-01, 7.85618966e-01],  # Schott N-PK51
    [5.16800155e-03, 1.66658798e-02, 1.38964129e02, 1.02960700e00, 1.88050600e-01, 7.36488165e-01],  # Schott N-PK52A
    [4.69824067e-03, 1.61818463e-02, 1.04374975e02, 8.87272110e-01, 4.89592425e-01, 1.04865296e00],  # Schott N-PSK3
    [7.06416337e-03, 2.33251345e-02, 9.74847345e01, 1.38121836e00, 1.96745645e-01, 8.86089205e-01],  # Schott N-PSK53A
    [1.19654879e-02, 5.90589722e-02, 1.35521676e02, 1.60865158e00, 2.37725916e-01, 1.51530653e00],  # Schott N-SF1
    [1.22241457e-02, 5.95736775e-02, 1.47468793e02, 1.62153902e00, 2.56287842e-01, 1.64447552e00],  # Schott N-SF10
    [1.31887070e-02, 6.23068142e-02, 1.55236290e02, 1.73759695e00, 3.13747346e-01, 1.89878101e00],  # Schott N-SF11
    [1.30512113e-02, 6.13691880e-02, 1.49517689e02, 1.69022361e00, 2.88870052e-01, 1.70451870e00],  # Schott N-SF14
    [1.16507014e-02, 5.97856897e-02, 1.32709339e02, 1.57055634e00, 2.18987094e-01, 1.50824017e00],  # Schott N-SF15
    [1.09019098e-02, 5.85683687e-02, 1.27404933e02, 1.47343127e00, 1.63681849e-01, 1.36920899e00],  # Schott N-SF2
    [1.26793450e-02, 6.02038419e-02, 1.45760496e02, 1.67780282e00, 2.82849893e-01, 1.63539276e00],  # Schott N-SF4
    [1.12547560e-02, 5.88995392e-02, 1.29141675e02, 1.52481889e00, 1.87085527e-01, 1.42729015e00],  # Schott N-SF5
    [1.41749518e-02, 6.40509927e-02, 1.77389795e02, 1.87543831e00, 3.73757490e-01, 2.30001797e00],  # Schott N-SF57
    [1.33714182e-02, 6.17533621e-02, 1.74017590e02, 1.77931763e00, 3.38149866e-01, 2.08734474e00],  # Schott N-SF6
    [1.47053225e-02, 6.92998276e-02, 1.61817601e02, 2.02459760e00, 4.70187196e-01, 2.59970433e00],  # Schott N-SF66
    [1.33714182e-02, 6.17533621e-02, 1.74017590e02, 1.77931763e00, 3.38149866e-01, 2.08734474e00],  # Schott N-SF6HT
    [1.14338344e-02, 5.82725652e-02, 1.33241650e02, 1.55075812e00, 2.09816918e-01, 1.46205491e00],  # Schott N-SF8
    [6.80282081e-03, 2.19737205e-02, 1.01513232e02, 1.17963631e00, 2.29817295e-01, 9.35789652e-01],  # Schott N-SK11
    [4.61716525e-03, 1.68859270e-02, 1.03736265e02, 9.36155374e-01, 5.94052018e-01, 1.04374583e00],  # Schott N-SK14
    [7.04687339e-03, 2.29005000e-02, 9.27508526e01, 1.34317774e00, 2.41144399e-01, 9.94317969e-01],  # Schott N-SK16
    [7.27191640e-03, 2.42823527e-02, 1.10377773e02, 1.28189012e00, 2.57738258e-01, 9.68186040e-01],  # Schott N-SK2
    [7.16874107e-03, 2.46455892e-02, 1.00886364e02, 1.32993741e00, 2.28542996e-01, 9.88465211e-01],  # Schott N-SK4
    [5.22730467e-03, 1.72733646e-02, 9.83594579e01, 9.91463823e-01, 4.95982121e-01, 9.87393925e-01],  # Schott N-SK5
    [8.23982975e-03, 3.33736841e-02, 1.06870822e02, 1.43060270e00, 1.53150554e-01, 1.01390904e00],  # Schott N-SSK2
    [9.20284626e-03, 4.23530072e-02, 1.06927374e02, 1.59222659e00, 1.03520774e-01, 1.05174016e00],  # Schott N-SSK5
    [8.69310149e-03, 4.21566593e-02, 1.11300666e02, 1.44857867e00, 1.17965926e-01, 1.06937528e00],  # Schott N-SSK8
    [6.76601657e-03, 2.30642817e-02, 8.90498778e01, 1.07715032e00, 1.68079109e-01, 8.51889892e-01],  # Schott N-ZK7
    [6.76601657e-03, 2.30642817e-02, 8.90498778e01, 1.07509891e00, 1.68895044e-01, 8.60503983e-01],  # Schott N-ZK7A
    [7.22141956e-03, 2.68216805e-02, 1.01702362e02, 1.18318503e00, 8.71756426e-02, 1.03133701e00],  # Schott P-BK7
    [9.38006396e-03, 3.60537464e-02, 8.64324693e01, 1.76003244e00, 2.48286745e-01, 1.15935122e00],  # Schott P-LAF37
    [7.15959695e-03, 2.33637446e-02, 8.83284426e01, 1.39324260e00, 4.18882766e-01, 1.04380700e00],  # Schott P-LAK35
    [1.00328203e-02, 3.87095168e-02, 9.45421507e01, 1.85543101e00, 3.15854649e-01, 1.28561839e00],  # Schott P-LASF47
    [9.99234757e-03, 3.87437988e-02, 9.58967681e01, 1.84910553e00, 3.29828674e-01, 1.30400901e00],  # Schott P-LASF50
    [9.88495571e-03, 3.78097402e-02, 9.78415430e01, 1.84568806e00, 3.39001600e-01, 1.32418921e00],  # Schott P-LASF51
    [1.68838419e-02, 7.16086325e-02, 1.18707479e02, 2.33300670e00, 4.52961396e-01, 1.25172339e00],  # Schott P-SF68
    [1.21696677e-02, 6.00710405e-02, 1.45651908e02, 1.62594647e00, 2.35927609e-01, 1.67434623e00],  # Schott P-SF69
    [1.16582670e-02, 5.82087757e-02, 1.30748028e02, 1.55370411e00, 2.06332561e-01, 1.39708831e00],  # Schott P-SF8
    [7.40877235e-03, 2.54563489e-02, 1.07751087e02, 1.31053414e00, 1.69376189e-01, 1.10987714e00],  # Schott P-SK57
    [7.36408831e-03, 2.55786047e-02, 1.06726060e02, 1.30536483e00, 1.71434328e-01, 1.10117219e00],  # Schott P-SK57Q1
    [7.20717498e-03, 2.45659595e-02, 1.02739728e02, 1.31678410e00, 1.71154756e-01, 1.12501473e00],  # Schott P-SK58A
    [7.84382378e-03, 2.87769365e-02, 1.05373397e02, 1.40790442e00, 1.43381417e-01, 1.16513947e00],  # Schott P-SK60
    [1.21481001e-02, 5.34549042e-02, 1.12174809e02, 1.55912923e00, 2.84246288e-01, 9.68842926e-01],  # Schott SF1
    [1.27534559e-02, 5.81983954e-02, 1.16607680e02, 1.61625977e00, 2.59229334e-01, 1.07762317e00],  # Schott SF10
    [1.36068604e-02, 6.15960463e-02, 1.21922711e02, 1.73848403e00, 3.11168974e-01, 1.17490871e00],  # Schott SF11
    [1.05795466e-02, 4.93226978e-02, 1.12405955e02, 1.40301821e00, 2.31767504e-01, 9.39056586e-01],  # Schott SF2
    [1.25502104e-02, 5.44559822e-02, 1.17652222e02, 1.61957826e00, 3.39493189e-01, 1.02566931e00],  # Schott SF4
    [1.11826126e-02, 5.08594669e-02, 1.12041888e02, 1.46141885e00, 2.47713019e-01, 9.49995832e-01],  # Schott SF5
    [1.33874699e-02, 5.79561608e-02, 1.21616024e02, 1.70579259e00, 3.44223052e-01, 1.09601828e00],  # Schott SF56A
    [1.43704198e-02, 5.92801172e-02, 1.21419942e02, 1.81651371e00, 4.28893641e-01, 1.07186278e00],  # Schott SF57
]

ALL_GLASS_NAMES = np.array(
    [
        "SiO2",
        "GeO2",
        "9.1% P2O2",
        "13.3% B2O3",
        "1.0% F",
        "16.9% Na2O : 32.5% B2O3",
        "ABCY",
        "HBL",
        "ZBG",
        "ZBLA",
        "ZBLAN",
        "5.2% B2O3",
        "10.5% P2O2",
        "N-BK7",
        "fused silica",
        "sapphire (ordinary)",
        "sapphire (extraordinary)",
        "MgF2 (ordinary)",
        "MgF2 (extraordinary)",
        "CaF2",
        "F2",
        "F5",
        "FK5HTi",
        "K10",
        "K7",
        "LAFN7",
        "LASF35",
        "LF5",
        "LF5HTi",
        "LLF1",
        "LLF1HT",
        "N-BAF1",
        "N-BAF4",
        "N-BAF5",
        "N-BAF5",
        "N-BAK1",
        "N-BAK2",
        "N-BAK4",
        "N-BALF",
        "N-BALF",
        "N-BASF",
        "N-BASF",
        "N-BK10",
        "N-F2",
        "N-FK5",
        "N-FK51",
        "N-FK58",
        "N-K5",
        "N-KF9",
        "N-KZFS",
        "N-KZFS",
        "N-KZFS",
        "N-KZFS",
        "N-KZFS",
        "N-LAF2",
        "N-LAF2",
        "N-LAF3",
        "N-LAF3",
        "N-LAF3",
        "N-LAF7",
        "N-LAK1",
        "N-LAK1",
        "N-LAK1",
        "N-LAK2",
        "N-LAK2",
        "N-LAK3",
        "N-LAK3",
        "N-LAK7",
        "N-LAK8",
        "N-LAK9",
        "N-LASF",
        "N-LASF",
        "N-LASF",
        "N-LASF",
        "N-LASF",
        "N-LASF",
        "N-LASF",
        "N-LASF",
        "N-LASF",
        "N-PK51",
        "N-PK52",
        "N-PSK3",
        "N-PSK5",
        "N-SF1",
        "N-SF10",
        "N-SF11",
        "N-SF14",
        "N-SF15",
        "N-SF2",
        "N-SF4",
        "N-SF5",
        "N-SF57",
        "N-SF6",
        "N-SF66",
        "N-SF6H",
        "N-SF8",
        "N-SK11",
        "N-SK14",
        "N-SK16",
        "N-SK2",
        "N-SK4",
        "N-SK5",
        "N-SSK2",
        "N-SSK5",
        "N-SSK8",
        "N-ZK7",
        "N-ZK7A",
        "P-BK7",
        "P-LAF3",
        "P-LAK3",
        "P-LASF",
        "P-LASF",
        "P-LASF",
        "P-SF68",
        "P-SF69",
        "P-SF8",
        "P-SK57",
        "P-SK57",
        "P-SK58",
        "P-SK60",
        "SF1",
        "SF10",
        "SF11",
        "SF2",
        "SF4",
        "SF5",
        "SF56A",
        "SF57",
    ]
)


def glass(number):
    """
    Return an array of Sellmeier coefficients for glass.

    The glasses all have a number.  This number is used to look up
    the Sellmeier coefficients for glass.  The number can be looked
    up with `ofiber.glass_name`.

    Use like this::
        num = ofiber.glass_index("SiO2")
        glass = ofiber.refraction.glass(num)       # SiO2
        n = ofiber.refraction.n(glass,632.8e-9)

    Args:
        number: glass number
    Returns:
        array of Sellmeier coefficients
    """
    return _GLASS[number]


def glass_name(number):
    """
    Look up the name for a specific glass number.

    A list of all possible names is in the array
    ofiber.refraction.ALL_GLASS_NAMES

    Args:
        number: glass number
    Returns:
        name of glass
    """
    return ALL_GLASS_NAMES[number]


def find_glass(name):
    """
    Look up the index of the glass with a particular name.

    The index of the first glass that matches the string is returned.
    Matching is case insensitive

    Args:
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
    Calculate Sellmeier coefficients for SiO_2 doped with GeO_2.

    The idea is that the glass a combination of silicon dioxide or
    germanium dioxide where x is the molar fraction of GeO_2. Thus
    (1-x) is the molar fraction of SiO_2.  The overall composition is
    x * GeO_2 : (1 - x) * SiO_2.

    Args:
        x: fraction of GeO_2 (0<=x<=1)

    Returns:
        Sellmeier coefficients for doped glass (array of six values)
    """
    SA = np.array([0.6961663, 0.4079426, 0.8974794])
    SL = np.array([0.0684043, 0.1162414, 9.896161])
    GA = np.array([0.80686642, 0.71815848, 0.85416831])
    GL = np.array([0.068972606, 0.15396605, 11.841931])
    a = (SL + x * (GL - SL)) ** 2
    b = abs(SA + x * (GA - SA))
    return np.concatenate([a, b])


def doped_glass_name(x):
    """
    Create a string the name describing the GeO_2 doped glass.

    Args:
        x: molar fraction of GeO_2 in the system
    Returns:
        string describing the doped glass
    """
    if x == 0:
        return r"SiO$_2$"
    if x == 1:
        return r"GeO$_2$"

    return r"%.2f GeO$_2$ : %.2f SiO$_2$" % (x, 1 - x)


def _sellmeier(b, c, lambda0):
    """
    Calculate the index of refraction using the Sellmeier equation.

    This is intended as a private method.

    Args:
        b: array of three Sellmeier Coefficients  [--]
        c: array of three Sellmeier Coefficients  [microns**2]
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
    Calculate the first derivative (wrt wavelength) of the Sellmeier equation.

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
        dy -= b[i] * c[i] / (lam2 - c[i]) ** 2  # 1/um**2

    dy *= lam / n1  # 1/um
    dy *= 1e6  # 1/m

    return dy


def _d2_sellmeier(b, c, lambda0):
    """
    Calculate the second derivative (wrt wavelength) of the Sellmeier equation.

    This is a private method.

    Args:
        b : array of three Sellmeier Coefficients  [--]
        c : array of three Sellmeier Coefficients  [microns**2]
        lambda0 : wavelength in vacuum             [m]

    Returns:
        returns the second derivative of the refractive index at lambda0 [1/m]
    """
    nn = _sellmeier(b, c, lambda0)  # index of refraction
    lam = lambda0 * 1e6  # needed because Sellmeier uses [um]
    lam2 = lam**2  # [um**2]

    dy = 0
    d2y = 0
    for i in range(3):
        dy = b[i] * c[i] / (lam2 - c[i]) ** 2  # 1/um
        d2y += b[i] * c[i] * (3 * lam2 + c[i]) / (lam2 - c[i]) ** 3  # 1/um**2

    total = d2y / nn - lam2 * dy**2 / nn**3  # 1/um**2
    total *= 1e12  # 1/m**2

    return total


def n(glass_coef, lambda0):
    """
    Calculate index of refraction for a glass at a wavelength.

    Args:
        glass_coef: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]

    Returns:
        index of refraction [--]
    """
    return _sellmeier(glass_coef[3:6], glass_coef[0:3], lambda0)


def dn(glass_coef, lambda0):
    """
    Calculate the first derivative of the refractive index w.r.t. wavelength.

    Args:
        glass_coef: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]

    Returns:
        the first derivative of index of refraction [1/m]
    """
    return _d_sellmeier(glass_coef[3:6], glass_coef[0:3], lambda0)


def d2n(glass_coef, lambda0):
    """
    Calculate the second derivative of the refractive index w.r.t. wavelength.

    Args:
        glass_coef: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]

    Returns:
        the second derivative of index of refraction [1/m**2]
    """
    return _d2_sellmeier(glass_coef[3:6], glass_coef[0:3], lambda0)


def n_group(glass_coef, lambda0):
    """
    Calculate group index of refraction at a wavelength.

    Args:
        glass_coef: array of Sellmeier coefficients obtained from glass(i)
        lambda0: wavelength in vacuum [m]

    Returns:
        group index of refraction [--]
    """
    return n(glass_coef, lambda0) - lambda0 * dn(glass_coef, lambda0)


def n_air(lambda0, temperature=15):
    """
    Calculate refractive index of air at atmospheric pressure.

    This follows the equation on page 4 of Smith, Modern Optical Engineering.

    Args:
        lambda0: wavelength in vacuum [m]
        temperature: degrees celsius
    Returns:
        index of refraction [--]
    """
    nu = 1 / (lambda0 * 1e6)
    n15 = 1e-8 * (8342.1 + 2406030 / (130 - nu**2) + 15996 / (38.9 - nu**2))
    if temperature == 15:
        return 1 + n15

    return 1 + 1.0549 * n15 / (1 + 0.00366 * temperature)


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
# extract glasses that are not SiO₂:GeO₂ mixtures
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
