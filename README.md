# CellCount_Nishijima_2024
# Fecal microbial load is a major determinant of gut microbiome variation and a confounder for disease associations

Suguru Nishijima, Evelina Stankevic, Oliver Aasmets, Thomas S. B. Schmidt, Naoyoshi Nagata, Marisa Isabell Keller, Pamela Ferretti, Helene Bæk Juel, Anthony Fullam, Shahriyar Mahdi Robbani, Christian Schudoma, Johanne Kragh Hansen, Louise Aas Holm, Mads Israelsen, Robert Schierwagen, Nikolaj Torp, Manimozhiyan Arumugam, Flemming Bendtsen, Charlotte Brøns, Cilius Esmann Fonvig, Jens-Christian Holm, Trine Nielsen, Julie Steen Pedersen, Maja Sofie Thiele, Jonel Trebicka, Elin Org, Aleksander Krag, Torben Hansen, Michael Kuhn, and Peer Bork, on behalf of the GALAXY and MicrobLiver Consortia

# R codes used to construct prediction models
This repository contains R codes to construct prediction models for fecal microbial load (i.e. total microbial cells / g, cell density) and evaluate their accuracy. The code uses species-level taxonomic profiles (mOTUs v2.5) and fecal microbial load of the GALAXY/MicrobLiver and MetaCardis study populations, which are contained in the `data` folder. Additionally, the repository contains all the R codes to generate figures with data that can be made public, except for the restricted Japense 4D and Estonian Microbiome datasets. 

To generate figures, it is first necessary to construct prediction models using construct_models.R and . 
