from lib.pcoa.pcoa_polar import calc_pcoa

def exe_pcoa():
    root = input("保存先のrootを入力") # root = "/home/nomu/work/Xenium/3_calc_min/result/min_vector/breast/20240702140629/2400_3400_2400_4000/"
                                    # mesh_calc_min.pyの返り値において、path_nameにあたる
    standard_gene = input("基準遺伝子を入力")  #  CDH1 など, # cdh1_file = root + standard_gene + "_dictionary.json" で使用
    coordinates = calc_pcoa(root, standard_gene)

    return coordinates
