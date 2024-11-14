# CellCount_Nishijima_2024
# Fecal microbial load is a major determinant of gut microbiome variation and a confounder for disease associations
# https://www.sciencedirect.com/science/article/pii/S0092867424012042
# https://doi.org/10.1016/j.cell.2024.10.022

Suguru Nishijima, Evelina Stankevic, Oliver Aasmets, Thomas S. B. Schmidt, Naoyoshi Nagata, Marisa Isabell Keller, Pamela Ferretti, Helene Bæk Juel, Anthony Fullam, Shahriyar Mahdi Robbani, Christian Schudoma, Johanne Kragh Hansen, Louise Aas Holm, Mads Israelsen, Robert Schierwagen, Nikolaj Torp, Manimozhiyan Arumugam, Flemming Bendtsen, Charlotte Brøns, Cilius Esmann Fonvig, Jens-Christian Holm, Trine Nielsen, Julie Steen Pedersen, Maja Sofie Thiele, Jonel Trebicka, Elin Org, Aleksander Krag, Torben Hansen, Michael Kuhn, and Peer Bork, on behalf of the GALAXY and MicrobLiver Consortia

# R codes used to construct prediction models
This repository contains R codes to construct prediction models for fecal microbial load (i.e. total microbial cells / g, cell density) and evaluate their accuracy. The code uses species-level taxonomic profiles (mOTUs v2.5) and fecal microbial load of the GALAXY/MicrobLiver (n = 1,894) and MetaCardis (n = 1,812) study populations, which are contained in the `data` folder. Additionally, the repository contains all the R codes to generate figures with data that can be made public, except for the restricted Japense 4D and Estonian Microbiome datasets. 

To generate figures, it is first necessary to construct prediction models and validate them using construct_models.R, internal_validation.R, external_validation_of_models.R. Since it takes hours to days to construct all the models, we also provide already constructed models in `out/model/` directroy. If you want to use those models, you can start with internal_validation.R and external_validation_of_models.R instead of using construct_models.R.
