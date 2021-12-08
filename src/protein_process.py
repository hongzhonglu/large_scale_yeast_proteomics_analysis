import math



def getSurfaceRatio(volume_ratio = 0.005/100, Vcell = 82, Scell=91.27):
    """
    The function is used to calculate the surface area ratio between organelle and cell.
    :param volume_ratio:
    :param Vcell:
    :param Scell:
    :return:
    """
    Vany = Vcell * volume_ratio
    Dany = 2 * (3 * Vany / (4 * math.pi)) ** (1 / 3)
    Sanym = 4 * math.pi * (Dany / 2) ** 2
    ratio = Sanym / Scell
    return ratio

def getGeneListFromLocation(gene_location_annotation, location):
    """
    The function is to extract gene list based on its compartment information
    :param gene_location_annotation:
    :param location:
    :return:
    """
    #location = 'mitochondrial envelope'
    gene_subset = gene_location_annotation[gene_location_annotation["GO_Name"]==location]
    gene_list = list(set(gene_subset["Systematic_name"].tolist()))
    return gene_list








