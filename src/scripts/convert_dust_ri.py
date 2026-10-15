#!/usr/bin/env python3.11
import os

# List of sites
sites = ["Algeria", "Arizona", "Atacama", "Australia",
         "Bodele", "Ethiopia", "Gobi", "Kuwait", 
         "Libya", "Mali", "Mauritania", "Morocco",
         "Namib-1", "Namib-2", "Niger", "Patagonia",
         "SaudiArabia", "Taklimakan", "Tunisia"]

for site in sites:
    fn = f"../geosmie/data/dust_composition/{site}_cri_selected.txt"
    if os.path.exists(fn):
        print(f"TRUE:  {fn}")
    else:
        print(f"FALSE: {fn}")

    fo = f"../geosmie/data/dust_composition/ri-du_{site}.wsv"

    with open(fn,"r") as f:
        data = f.readlines()
    f.close

    with open(fo,"w") as f:
        f.write("# lambda[um] m_real m_imaginary\n")
        i = 0
        for line in data:
            if(i<2):
                i += 1
                continue
            i += 1
            items = line.split()
            f.write(f"{items[0]} {items[1]} -{items[2]}\n")
    f.close()
