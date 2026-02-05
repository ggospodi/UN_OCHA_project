# 🇳🇵 Nepal Earthquake Analysis: Aid Relief Coverage, Severity Index & Aid Capacity Index

**A data science collaboration with the United Nations Office for the Coordination of Humanitarian Affairs (UN OCHA)**

*September–December 2015*

---

## Overview

This project analyzes the humanitarian aid response to the **April 2015 Nepal earthquake** (magnitude 7.8), which killed nearly 9,000 people, injured over 22,000, and displaced millions. Working directly with UN OCHA and the Humanitarian Data Exchange (HDX), our international team of data scientists built network-based models to evaluate how effectively aid reached the ~800 Village Development Committee (VDC) districts across Nepal.

The core contributions include:

- **Agency-VDC Aid Network** — a bipartite network mapping ~80 aid agencies to ~800 VDC districts, revealing the structure and gaps in aid delivery
- **Severity-Based Need Index** — a composite model estimating relative need for aid using hazard exposure, vulnerability, distance to epicenter, and healthcare facility density
- **Aid Capacity Index** — a novel metric capturing each region's potential to receive adequate aid coverage, weighted by agency hub scores in the network

> **Key finding:** Aid distribution need is largely determined by pre-existing conditions (vulnerability, infrastructure, healthcare access) and can therefore be anticipated and prepared for *before* disaster strikes.

## Presentation

The full analysis deck (17 slides) is available on SlideShare:

🔗 **[Nepal Disaster Analysis — SlideShare](https://www.slideshare.net/slideshow/nepal-disaster-analysis/61321723)**

## Methodology

### Need for Aid Model

We define a **Severity Index** following the [INFORM Index](http://www.inform-index.org/InDepth/Methodology) framework:

```
Severity = (Hazard × Exposure × Vulnerability) ^ (1/3)
```

Need for aid at each VDC is then modeled as:

```
Need ∝ Severity / (Distance_to_Epicenter^(1/3) × Healthcare_Facilities^(1/3))
```

This produces a relative ranking across all VDCs:
- **Low Need:** 64.1% of VDCs
- **Medium Need:** 33.0% of VDCs
- **High Need:** 2.9% of VDCs

### Agency-VDC Aid Network

The bipartite network connects aid agencies to VDC districts they served. Key statistics:

| Metric | Min | Median | Mean | Max |
|--------|-----|--------|------|-----|
| Aid instances per agency | 1 | 18 | 40.6 | 256 |
| VDCs served per agency | 1 | 9 | 17.7 | 127 |
| Agencies per VDC (shared) | 0 | 11 | 13 | 53 |

### Aid Coverage Findings

- Some VDCs received aid from up to **12 agencies**, while the average was just **2**
- Only ~**70%** of Medium-to-High Need VDCs received any aid
- Only ~**40%** of Low Need VDCs received aid per available records
- Significant gaps exist where high-need areas had single-agency coverage or no coverage at all

### Aid Capacity Index

```
Aid Capacity Index = f(Σ VDC_Degree × Agency_Hub_Score)
```

A weighted function of each VDC's degree in the aid network, scaled by the hub authority score of contributing agencies. This metric captures the *potential* for appropriate aid coverage should a disaster occur, enabling proactive resource positioning.

## Data Sources

| Source | Description |
|--------|-------------|
| [Humanitarian Data Exchange (HDX)](https://data.humdata.org/) | 600+ indicators from 30+ sources — the core humanitarian dataset |
| [INFORM Index](http://www.inform-index.org) | Hazard, exposure, and vulnerability indices |
| [WorldPop](http://www.worldpop.org.uk) | Population density and distribution data |
| [OSGeo / QGIS](http://www.osgeo.org) | Geospatial mapping of VDC boundaries |
| [ICIMOD](http://www.icimod.org/nepalearthquake2015) | Nepal earthquake impact assessment |

## Repository Structure

```
├── data/                          # Source datasets (HDX, INFORM, WorldPop)
├── 24hearthquakes.R               # Earthquake event analysis
├── Aid_And_Need.R                 # Aid vs. need comparison
├── Aid_And_Severity_table.R       # Severity tabulations
├── Aid_Severity_Model.R           # Primary severity model
├── Aid_Severity_Model2.R          # Model iteration
├── Aid_Severity_Model_Tables.R    # Model output tables
├── Disaster_Aid_Network.R         # Agency-VDC network construction
├── Disaster_Aid_Network.Rmd       # R Markdown report — aid network
├── Disaster_Aid_Network.html      # Rendered HTML report
├── Disaster_Aid_Network_Projections.R
├── Disaster_Aid_Network_Projections_hlcit.R
├── Disaster_Aid_Network_test.R
├── Disaster_Aid_Severity_Analysis.R
├── Displacement_Network.R         # Population displacement network
├── Displacement_Network.Rmd       # R Markdown report — displacement
├── Displacement_Network.html      # Rendered HTML report
├── Displacement_Severity_Analysis.R
├── Nepal_Network_Analysis_bkp.R   # Backup of network analysis
├── UN_OCHA_project.Rproj          # RStudio project file
└── README.md
```

> **Note:** This analysis was conducted in R (2015). The code is preserved as a reference implementation. The rendered HTML reports (`Disaster_Aid_Network.html`, `Displacement_Network.html`) contain the full interactive visualizations.

## Key Conclusions

1. **Humanitarian data initiatives** like HDX enable impactful, data-driven disaster response at global scale
2. **Aid distribution gaps** are significant — coordination across agencies is critical to avoid duplication in some areas and neglect in others
3. **Pre-existing conditions** (vulnerability, infrastructure, healthcare access) are the strongest predictors of aid need, meaning preparedness planning can be done *before* disasters strike
4. **Aid coverage follows a non-linear pattern** — simple geographic proximity to epicenter is insufficient for resource allocation
5. **The Aid Capacity Index** provides a framework for anticipating regional aid absorption capacity, enabling proactive positioning of resources

## Team

| Name | Role | Affiliation |
|------|------|-------------|
| **Georgi D. Gospodinov, PhD** | Lead Data Scientist | Walmart Stores Inc. |
| **Chris Shannon** | Sr. Data Scientist | Tata Consultancy Services |
| **Aidan McGuire** | HDX Product Manager | UN OCHA / ScraperWiki |
| **Javier Teran** | Statistician | United Nations OCHA |
| **Andrej Verity** | Information Management Officer | United Nations OCHA |

*Plus an international group of data science volunteers.*

## Context

This project was part of the broader [Humanitarian Data Exchange](https://data.humdata.org/) initiative by UN OCHA, which aims to make humanitarian data easy to find and use for analysis. The Nepal earthquake response involved a $415 million initial UN appeal, with the largest allocations going to food ($128M), health ($75M), and shelter ($50M) for the first three months of operations.

## License

This project was conducted as a volunteer humanitarian data science effort in collaboration with UN OCHA. The analysis and code are shared for educational and research purposes.

---

*Analysis conducted September–December 2015. Repository preserved as a reference implementation.*
