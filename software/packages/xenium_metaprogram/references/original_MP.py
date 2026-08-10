import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

class MalignantMetaProgrammer:
    """
    Advanced pipeline for discovering conserved meta-programs.
    Now includes Patient-Centering to suppress clonal signals.
    """
    
    def __init__(self, adata, sample_col='patient_id', layer='counts'):
        self.adata = adata
        self.sample_col = sample_col
        self.layer = layer
        
        # Internal state
        self.all_programs = [] 
        self.similarity_matrix = None
        self.meta_programs = {}
        self.program_assignments = None

    def preprocess_centered_pos(self, new_layer='centered_pos'):
        """
        [Step 0] Patient-Centering: V_centered = max(0, V - mean(V_patient))
        Crucial for ESCC/Malignant cells to remove the CNV/Clonal baseline.
        """
        print(f"✨ Neutralizing patient clonal signals in layer: {self.layer}...")
        
        # Extract data from specified layer
        v_matrix = self.adata.layers[self.layer].copy()
        if hasattr(v_matrix, 'toarray'): v_matrix = v_matrix.toarray()
        
        # Center per patient
        for sample in self.adata.obs[self.sample_col].unique():
            idx = self.adata.obs[self.sample_col] == sample
            # Subtract the average expression profile of this patient
            v_matrix[idx] -= v_matrix[idx].mean(axis=0)
        
        # Threshold at zero (NMF requirement)
        v_matrix[v_matrix < 0] = 0
        
        # Store and update class state
        self.adata.layers[new_layer] = v_matrix
        self.layer = new_layer
        print(f"✅ Created '{new_layer}'. Data is now relative to patient baseline.")

    def run_pipeline(self, k_range=(7, 14), top_n_genes=30, n_mps=15, min_patients=3, 
                     center_first=True, plot_results=True):
        """Orchestrates: Centering -> Multi-K NMF -> Clustering -> Scoring."""
        
        # 🔥 CRITICAL FIX: Clear previous runs so we don't accumulate thousands of programs
        self.all_programs = [] 
        self.meta_programs = {}
        
        if center_first:
            self.preprocess_centered_pos()

        print(f"\n[Step 1/4] Running Multi-K NMF (k={k_range[0]} to {k_range[1]})...")
        for k in range(k_range[0], k_range[1] + 1):
            print(f"  > Decomposing at rank k={k}...")
            self._run_sample_nmf(k=k, top_n_genes=top_n_genes)
        
        print(f"[Step 2/4] Building Similarity Matrix for {len(self.all_programs)} programs...")
        self._compute_similarity()
        
        print(f"[Step 3/4] Clustering Meta-Programs...")
        self._identify_meta_programs(n_clusters=n_mps, min_patient_threshold=min_patients)
        
        print(f"[Step 4/4] Scoring cells and finalizing...")
        self.score_cells()

        self._print_hallmark_report()
        
        print(f"✅ Pipeline Complete: {len(self.meta_programs)} MPs identified.")
        if plot_results:
            self.plot_nature_style_heatmap()

    def _get_nmf_loadings(self, sample_id, k):
        sample_adata = self.adata[self.adata.obs[self.sample_col] == sample_id].copy()
        data = sample_adata.layers[self.layer]
        if hasattr(data, 'toarray'): data = data.toarray()
        
        # Use nndsvd for better initialization with sparse/centered data
        model = NMF(n_components=k, init='nndsvd', random_state=42, max_iter=1000)
        model.fit(data)
        return model.components_, sample_adata.var_names

    def _run_sample_nmf(self, k, top_n_genes):
        samples = self.adata.obs[self.sample_col].unique()
        for sample_id in samples:
            try:
                H, feature_names = self._get_nmf_loadings(sample_id, k)
                for i in range(k):
                    top_indices = H[i].argsort()[-top_n_genes:][::-1]
                    self.all_programs.append({
                        'sample': sample_id,
                        'genes': set(feature_names[top_indices]),
                        'program_id': f"{sample_id}_k{k}_P{i}"
                    })
            except Exception as e:
                print(f"  ⚠️ Skipping {sample_id} at k={k} due to error (likely too few cells).")

    def _compute_similarity(self):
        n = len(self.all_programs)
        sim_mat = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                s1, s2 = self.all_programs[i]['genes'], self.all_programs[j]['genes']
                intersect = len(s1.intersection(s2))
                union = len(s1.union(s2))
                sim_mat[i, j] = sim_mat[j, i] = intersect / union if union > 0 else 0
        self.similarity_matrix = pd.DataFrame(sim_mat, 
                                              index=[p['program_id'] for p in self.all_programs],
                                              columns=[p['program_id'] for p in self.all_programs])

    def _identify_meta_programs(self, n_clusters, min_patient_threshold):
        Z = linkage(squareform(1 - self.similarity_matrix), method='average')
        self.program_assignments = fcluster(Z, t=n_clusters, criterion='maxclust')
        prog_df = pd.DataFrame(self.all_programs)
        prog_df['cluster'] = self.program_assignments
        
        for clust_id in range(1, n_clusters + 1):
            subset = prog_df[prog_df['cluster'] == clust_id]
            if subset['sample'].nunique() < min_patient_threshold: continue
            
            all_genes = [gene for gset in subset['genes'] for gene in gset]
            counts = pd.Series(all_genes).value_counts()
            # Consensus: gene must appear in > 20% of sample-programs in this cluster
            consensus = counts[counts >= len(subset) * 0.30].index.tolist()
            self.meta_programs[f"MP_{clust_id}"] = consensus

    def score_cells(self):
        for mp_name, genes in self.meta_programs.items():
            sc.tl.score_genes(self.adata, gene_list=genes, score_name=mp_name)

    def _print_hallmark_report(self, markers=['MKI67', 'TOP2A', 'VEGFA', 'MMP14', "COL17A1", "LAMA3"]):
        print("\n--- 🧬 Hallmark Gene Report ---")
        for marker in markers:
            found_in = [mp for mp, genes in self.meta_programs.items() if marker in genes]
            status = f"✅ found in {found_in}" if found_in else "❌ Not found"
            print(f"{marker:8}: {status}")

    def plot_nature_style_heatmap(self, vmax=0.3):
        """Visualizes the similarity matrix with strict index alignment."""
        if self.similarity_matrix is None: return

        # 1. Create order dataframe and ensure it matches the matrix indices exactly
        order_df = pd.DataFrame({
            'program_id': self.similarity_matrix.index,
            'cluster': self.program_assignments
        }).sort_values('cluster').reset_index(drop=True)
        
        # 2. Reorder the matrix
        reordered_mat = self.similarity_matrix.loc[order_df['program_id'], order_df['program_id']]
        
        # 3. Create cluster color mapping
        unique_clusters = sorted(order_df['cluster'].unique())
        palette = sns.color_palette("husl", len(unique_clusters))
        cluster_color_map = dict(zip(unique_clusters, palette))
        
        # 🔥 FIX: Convert to LIST to avoid Seaborn indexing errors
        network_colors = [cluster_color_map[c] for c in order_df['cluster']]

        # 4. Use Clustermap
        g = sns.clustermap(
            reordered_mat,
            row_cluster=False, 
            col_cluster=False,
            row_colors=network_colors, 
            col_colors=network_colors,
            cmap="mako",
            vmax=vmax,
            figsize=(10, 10),
            xticklabels=False,
            yticklabels=False
        )
        plt.show()

    def test_sample_run(self, sample_id=None, k=8, top_n_genes=30):
        if sample_id is None: sample_id = self.adata.obs[self.sample_col].unique()[0]
        H, names = self._get_nmf_loadings(sample_id, k)
        plt.figure(figsize=(12, 6))
        for i in range(min(k, 4)): # Plot first 4 for brevity
            top = H[i].argsort()[-top_n_genes:][::-1]
            plt.subplot(2, 2, i+1)
            plt.bar(names[top], H[i][top], color='gray')
            plt.xticks(rotation=90, fontsize=7); plt.title(f"Program {i}")
        plt.tight_layout(); plt.show()

    def plot_spatial_samples(self, mp_ids=None, spot_size=30, cmap='magma'):
        """
        Generates a master grid of spatial plots. 
        Rows = Samples, Columns = Meta-Programs.
        """
        # 1. Identify which MPs to plot
        mps_to_plot = mp_ids if mp_ids else list(self.meta_programs.keys())
        samples = self.adata.obs[self.sample_col].unique()
        
        n_samples = len(samples)
        n_mps = len(mps_to_plot)
        
        # 2. Setup the figure
        fig, axes = plt.subplots(n_samples, n_mps, figsize=(n_mps * 5, n_samples * 5))
        
        # Handle the case of a single sample or single MP (matplotlib axis array shape)
        if n_samples == 1 and n_mps == 1: axes = np.array([[axes]])
        elif n_samples == 1: axes = axes[np.newaxis, :]
        elif n_mps == 1: axes = axes[:, np.newaxis]

        print(f"🎨 Generating spatial grid for {n_samples} samples and {n_mps} programs...")

        for i, s_id in enumerate(samples):
            # Subset adata to just this patient
            sample_adata = self.adata[self.adata.obs[self.sample_col] == s_id]
            
            for j, mp_id in enumerate(mps_to_plot):
                # We use sc.pl.spatial if the coordinates are in obsm['spatial']
                sc.pl.spatial(
                    sample_adata, 
                    color=mp_id, 
                    cmap=cmap, 
                    spot_size=spot_size,
                    ax=axes[i, j], 
                    show=False, 
                    title=f"{s_id} | {mp_id}",
                    frameon=False
                )
                # Remove axis labels for a cleaner "Nature-style" look
                axes[i, j].set_xlabel("")
                axes[i, j].set_ylabel("")

        plt.tight_layout()
        plt.show()



#===========
# usage examples	

epi_adata =sc.read_h5ad("../../data/combined_adata/combined_all_patients_epi_analyzed_final.h5ad")
pre_treat_epi_adata = epi_adata[epi_adata.obs["sample_timepoint"]=="pre",]
pre_treat_epi_adata = pre_treat_epi_adata[pre_treat_epi_adata.obs["final_cell_type"]=="malignant"]
pre_treat_epi_adata

mp_analyzer = MalignantMetaProgrammer(pre_treat_epi_adata, sample_col="sampleID")

mp_analyzer.run_pipeline(
    k_range=(7,15),
    top_n_genes=50
    )