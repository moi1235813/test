from foamlib import FoamCase
import time


case1 = FoamCase("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/temp/T300")
case2 = FoamCase("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/temp/T500")
#case3 = FoamCase("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/temp/T1000")
#case4 = FoamCase("/home/jatant/OpenFOAM/jatant-14/run/jet2/boom/Burke/temp/T1300")

case_list = [case1, case2]
#case1.run(log=False)



"""for caseN in case_list:
    with caseN.control_dict as file:
        file["functions"]["HRR"]["writeInterval"] = 1e-07
        file["functions"]["specieFractions"]["writeInterval"] = 1e-07
        print(file["functions"]["HRR"]["writeInterval"])
        print(file["functions"]["specieFractions"]["writeInterval"])"""


for caseN in case_list:
    startTime = time.time()
    caseN.run()
    print(f"{caseN} completed in {time.time() - startTime} seconds")
    #with caseN["0"]["U"] as field:
        #print(field.boundary_field["inlet"])