# ValApp-for-singlecell-seurat-clustering
I built this app so my team members can review clustering and annotation at their own pace. It lets users explore different resolutions and view outputs like stacked bar plots, DEGs, violin, and feature plots. Users can also search for genes or upload gene sets for analysis.

## 🧬 Seurat Explorer: Interactive Single-Cell Analysis App

An interactive Shiny-based web application for exploring and analyzing single-cell RNA-seq (scRNA-seq) data using the Seurat framework. This app enables intuitive visualization, clustering, and downstream analysis without requiring extensive coding.

⸻

### 🚀 Features
	•	📊 Upload and explore Seurat objects
	•	🔍 Interactive UMAP / t-SNE visualization
	•	🧬 Cluster identification and annotation
	•	📈 Differential expression (DEG) analysis
	•	🌈 Customizable plots (FeaturePlot, Violin, DotPlot)
	•	📥 Export results (tables, plots, DEG lists)
	•	🔗 Integrated pipeline for streamlined workflows

⸻

### 🖥️ Live App

Access the app here:

👉 http://implk-01:3838/

⚠️ Note: Ensure you are on the same network or have proper access permissions.

⸻

### 📦 Installation

1. Clone the repository

git clone https://github.com/Valisha/ValApp-for-singlecell-seurat-clustering.git
cd ValApp-for-singlecell-seurat-clustering

2. Install dependencies

Make sure you have R (>= 4.4) installed.

install.packages(c("shiny", "Seurat", "dplyr", "ggplot2", "patchwork", "circlize", "openxlsx", "readxl", "DT", "tibble", "scales"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")


⸻

## ▶️ Running the App

Option 1: Using script

./start_app.sh

Option 2: Run manually in R

shiny::runApp()


⸻

## ⚙️ Deployment (Shiny Server)

To deploy on a Shiny Server (CentOS/RHEL):
	1.	Place the app folder in:

/var/www/shiny-server/

	2.	Restart Shiny Server:

sudo systemctl restart shiny-server

	3.	Access via:

http://<server-ip>:3838/<app-folder>/


⸻

## 📁 Project Structure

├── app.R / ui.R / server.R
├── start_app.sh
├── renv/                # Reproducible environment
├── data/                # Input datasets
├── www/                 # Static assets
└── README.md


⸻

## 📊 Data Input
	•	Accepts Seurat objects (.rds)
	•	Ensure metadata and clustering are precomputed for best performance

⸻

## 📈 Output
	•	Differential expression tables (downloadable)
	•	Publication-ready plots
	•	Cluster annotations

⸻

## 🛠️ Troubleshooting

Common Issues

1. Package installation errors
	•	Ensure system dependencies (e.g., libxml2, icu) are installed

2. App not loading
	•	Check logs:

sudo journalctl -u shiny-server

3. Upload size limits
	•	Modify Shiny Server config:

/etc/shiny-server/shiny-server.conf

Add:

app_init_timeout 300;


⸻

### 🤝 Contributing

Contributions are welcome!
	1.	Fork the repo
	2.	Create a new branch
	3.	Submit a pull request

⸻

### 📜 License

This project is licensed under the MIT License.

⸻

### 👤 Author

Valisha Shah

⸻

#### ⭐ Acknowledgements
	•	Seurat Team
	•	R Shiny Community
	•	Bioconductor Project

⸻

#### 💡 Future Enhancements
	•	Integration with TCR-seq analysis
	•	Advanced trajectory analysis
	•	Multi-omics (CITE-seq) support

⸻

If you find this useful, consider giving the repo a ⭐!
