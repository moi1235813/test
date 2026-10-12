import foamlib as fl


case1 = fl.FoamCase("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/velo/U10")

case2 = case1.clone("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/velo/U5")
case2 = case1.clone("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/velo/U15")
case2 = case1.clone("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/velo/U20")
case2 = case1.clone("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/velo/U30")