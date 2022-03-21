# /***
# * BLASTphylo
# * Normalisation of the NCBI taxonomy
# *
# * @author jennifer mueller
# *
# ***/

# used libraries
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import json
import os

# general variables
normalisation_ncbi_taxonomy_dic = {}
searchbar_entries = []
count_taxa = 0
count_species = 0
# bacteria_taxa = ncbi.get_descendant_taxa(2, intermediate_nodes=True)
# bacteria_ranks = ncbi.get_rank(bacteria_taxa)
# bacteria_tree = ncbi.get_topology(bacteria_taxa)

# bacteria
#species_tree = ncbi.get_descendant_taxa(2, intermediate_nodes=True)
# archaea
species_tree = ncbi.get_descendant_taxa(2157, intermediate_nodes=True)

species_tree = ncbi.get_topology(species_tree)

########################################################################################################################
#                                                                                            normalisation NCBI taxonomy
########################################################################################################################

''' Generate json file for the normalisation of the NCBI taxonomy 
    Function count the number of taxon with the lowest rank (strain or species) for each taxon in the NCBI taxonomy as 
    reference for the normalisation
'''
# def normalisation_ncbi_taxonomy(tree):
#     global count_taxa
#
#     if len(tree.get_children()) == 0:
#         normalisation_ncbi_taxonomy_dic[tree.sci_name] = [tree.name, 1]
#         searchbar_entries.append({'id': count_taxa, 'name': tree.sci_name + ' (txid: ' + tree.name + ')', 'menuInput': tree.name})
#         count_taxa += 1
#     else:
#         [normalisation_ncbi_taxonomy(child) for child in tree.get_children()]
#         normalisation_ncbi_taxonomy_dic[tree.sci_name] = [tree.name, sum([normalisation_ncbi_taxonomy_dic[child.sci_name][1] for child in tree.get_children()])]
#         print([tree.sci_name, normalisation_ncbi_taxonomy_dic[tree.sci_name]])
#         searchbar_entries.append({'id': count_taxa, 'name': tree.sci_name + ' (txid: ' + tree.name + ')', 'menuInput': tree.name})
#         count_taxa += 1


def normalisation_ncbi_taxonomy_species(tree):
    global count_species
    #print(ncbi.get_rank([tree.name])[int(tree.name)])
    if (ncbi.get_rank([tree.name])[int(tree.name)]=='genus'):
        normalisation_ncbi_taxonomy_dic[tree.sci_name] = [tree.name, 1]
        count_species += 1
    else:
        if ncbi.get_rank([tree.name])[int(tree.name)] not in ['strain', 'species group', 'species']:
            [normalisation_ncbi_taxonomy_species(child) for child in tree.get_children()]
            normalisation_ncbi_taxonomy_dic[tree.sci_name] = [tree.name, sum(
                [normalisation_ncbi_taxonomy_dic[child.sci_name][1] for child in tree.get_children()])]
        else:
            normalisation_ncbi_taxonomy_dic[tree.sci_name] = [tree.name, 0]


def main():
    print("Program: Generate NCBI Taxonomy normalisation")
    print("Author: jennifer mueller")

    actual_dir = os.getcwd()
    #normalisation_ncbi_taxonomy(bacteria_tree)

    # with open(actual_dir + '/src/data/ncbi_normalisation.json', 'w') as outfile:
    #     json.dump(normalisation_ncbi_taxonomy_dic, outfile, indent=4)
    #
    # with open(actual_dir + '/src/data/searchbar_entries.json', 'w') as outfile:
    #     json.dump(searchbar_entries, outfile, indent=4)
    # print('complete')
    normalisation_ncbi_taxonomy_species(species_tree)
    with open('/Users/zabel/projects/blastphylo/BLASTphylo/src/data/archaea_ncbi_normalisation_genus.json', 'w') as outfile:
        json.dump(normalisation_ncbi_taxonomy_dic, outfile, indent=4)


if __name__ == "__main__":
    main()