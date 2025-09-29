# SoN_2023_Mammals
This repository contains the essential scripts and data needed to recreate the analyses and figures for the mammals group in the State of Nature Report 2023 by Hamaarag.

Important note - The mammal data in the Inner Sand unit was analyzed but no significant trends were found, thus the analysis isn't included in the repository.

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This repository contains the **essential scripts and data** needed to **recreate the analyses and figures** for the **mammal group** in the [**State of Nature Report 2023**](https://hamaarag.org.il/report/%d7%93%d7%95%d7%97-%d7%9e%d7%a6%d7%91-%d7%94%d7%98%d7%91%d7%a2-2022-%d7%9b%d7%a8%d7%9a-%d7%9e%d7%92%d7%95%d7%95%d7%9f-%d7%91%d7%99%d7%95%d7%9c%d7%95%d7%92%d7%99/) by Hamaarag.

---

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## Contents

This repository is organized into the following main directories and files:  

- **`R/`** – Contains R script for data analysis and visualization and html copies of scripts for comparison. 

- **`data/`** – Includes input datasets required for the analyses.

- **`output/`** – Stores the output of the analyses. Each unit that was analyzed has it's own folder

- **Project Configuration**:
  - **`SoN_2023_Mammals.Rproj`** – RStudio project file.
  - **`renv/`** & **`renv.lock`** – Dependency management files for reproducibility.
  - **`.gitignore`** & **`.Rprofile`** – Configuration and environment settings.  

- **Documentation**:
  - **`README.md`** – This documentation file.
  - **`LICENSE`** – License details.  

---

## How to Use

To reproduce the analyses and figures, follow these steps:

### 1. Clone the repository
Run the following command in a **Bash terminal** (e.g., Git Bash, macOS/Linux terminal, or Windows Command Prompt with Git installed):
```sh
git clone https://github.com/Hamaarag/SoN_2023_Mammals.git
cd SoN_2023_Mammals
```

### 2. Open the project in RStudio
- Open **RStudio** and use the menu:  
  `File` → `Open Project...` → Select `SoN_2023_Mammals.Rproj`
  
Alternatively, you can open the project directly from the command line (Mac/Linux):
```sh
open SoN_2023_Mammals.Rproj
```
or on Windows (if RStudio is in the PATH):
```sh
start SoN_2023_Mammals.Rproj
```

### 3. Install required R packages
Inside **RStudio**, install the necessary dependencies using `renv`:
```r
install.packages("renv")  # Install renv if not already installed. If running into problems, install version 1.0.7 using: `install.packages("renv", version = "1.0.7")`
renv::restore()           # Restore package environment
```

### 4. Run analyses
To execute the analyses and generate figures, run the RMarkdown scripts in the **R/** directory.

#### **Option 1: Run interactively in RStudio**
Open any `.Rmd` file (e.g., `R/Analyzing the Reptile data for State of Nature Report 2023.Rmd`) and click **"Knit"** to execute the script and generate outputs.

#### **Option 2: Run from the command line**
To render a specific RMarkdown file programmatically:
```r
rmarkdown::render("R/Analyzing the Mammal data for State of Nature Report 2023.Rmd")
```

### 5. View the results
- Figures and tables (if generated) will also be saved in `output/`, under a designated folder (the folder will be created
  automatically).
- The knitted reports (`.html`) will appear in the same directory as the executed `.Rmd` file.

---

## Additonal notes regarding environment

 - R version:4.3.2  
 - RStudio version: 2025.05.1 
 - OS: Windows 10 x64  
 - The Rmd files are set to use the "Almoni ML v5 AAA" font. If you do not have this font installed, you may need to change the font in the Rmd files.  
 
---

## Citation

If you use this repository or its contents, please cite **both** the main report and this repository.

### **State of Nature Report 2023 – Biodiversity**
Grossbard S, Renan I (Eds). 2024. *State of Nature Report 2023 – Biodiversity*.  
Hamaarag – Israel’s National Ecosystem Assessment Program, Steinhardt Museum of Natural History, Tel Aviv University.  

### **SoN_2023_Mammals Repository**
Hamaarag. *State of Nature Report 2023 - Mammals*. Public Repository, 2024.  
Available at: [https://github.com/Hamaarag/SoN_2023_Reptiles.git](https://github.com/Hamaarag/SoN_2023_Reptiles.git).  

---

## Contact

For questions or feedback, contact **Hamaarag** at [info@hamaarag.org.il](mailto:info@hamaarag.org.il).

---

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
