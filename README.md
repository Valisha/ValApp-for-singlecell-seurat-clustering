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

⚠️ Note: Only works for salk employees 😅

⸻

### Github Pages

Access the app here:

👉 https://valisha.github.io/ValApp-for-singlecell-seurat-clustering/

⸻

### 📦 Installation

#### Clone the repository

```
git clone https://github.com/Valisha/ValApp-for-singlecell-seurat-clustering.git
cd ValApp-for-singlecell-seurat-clustering
```

After cloning, there are multiple ways to run the server locally on your computer. In each case, once the the app is running, it should be available at https://0.0.0.0:3838/ .

#### Option 1: Install and run with Pixi

Pixi is an environment manager that ensures that the same set of package versions is used across installations. If you don't have Pixi yet, it can be obtained here: https://pixi.prefix.dev/latest/installation/ 

Once you have Pixi, dependencies can be installed and this app can be started using:
```
pixi run start
```

#### Option 2: Install and run with Docker

Docker provides a container that holds both an operating system and the dependencies of this app. You can obtain Docker Desktop here: https://www.docker.com/products/docker-desktop/ 

Once installed, start the Docker app. Then, this app can be started using:
```
docker compose up -d
```
The first time the app is started, it may take a couple of minutes for dependencies to be installed within the container.

The app can be stopped using
```
docker compose down
```

#### Option 3: Install dependencies and run in your local R environment

Make sure you have R (>= 4.4) installed, then add additional dependencies with:

```
install.packages(c("shiny", "Seurat", "dplyr", "ggplot2", "patchwork", "circlize", "openxlsx", "readxl", "DT", "tibble", "scales"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

Once installed, the app can then be run locally using:
```
./start_app.sh
```
Or manually within R or Rstudio with:
```
shiny::runApp()
```

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
