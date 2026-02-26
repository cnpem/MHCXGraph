# utils/vis_tracer.py
from pathlib import Path
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

try:
    import imageio.v2 as imageio  # para GIF
except Exception:
    imageio = None

class TraversalTracer:
    """
    Coleta estados durante a busca e renderiza um GIF ou MP4 por âncora.
    Estados gravados por frame:
      current  : nó sendo processado no momento
      chosen   : tuple de nós já aceitos no ramo
      frontier : tuple dos nós na fronteira
      accepted_edges_last_leaf : edges do último subgrafo aceito neste âncora, p destacar
    """

    def __init__(self, out_dir, fmt="gif", fps=12, sample_every=25, max_frames=3000,
                 dpi=110, edge_alpha=0.45, seed=42, enabled=True):
        self.enabled = enabled
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.fmt = fmt  # "gif" ou "mp4"
        self.fps = fps
        self.sample_every = max(1, int(sample_every))
        self.max_frames = max_frames
        self.dpi = dpi
        self.edge_alpha = edge_alpha
        self.seed = seed

        self.reset()

    def reset(self):
        self.anchor = None
        self.frames = []
        self.step = 0
        self.pos = None
        self.G = None
        self.accepted_edges_last_leaf = set()

    def start_anchor(self, anchor, A_sub, nodes_map=None):
        """
        anchor: int id do âncora
        A_sub : matriz de adjacência binária do componente que está sendo explorado
        nodes_map: np.ndarray com ids verdadeiros se você estiver explorando um subgrafo; se None, assume range
        """
        if not self.enabled:
            return
        self.reset()
        self.anchor = int(anchor)

        # constrói grafo networkx uma vez
        n = A_sub.shape[0]
        if nodes_map is None:
            nodes_map = np.arange(n, dtype=int)
        G = nx.Graph()
        # adiciona vértices com seu id real
        for i in range(n):
            G.add_node(int(nodes_map[i]))
        # adiciona arestas
        ii, jj = np.where(np.triu(A_sub, 1))
        for i, j in zip(ii, jj):
            G.add_edge(int(nodes_map[i]), int(nodes_map[j]))
        self.G = G

        # layout fixo por âncora
        self.pos = nx.spring_layout(G, seed=self.seed, dim=2)

    def tick(self, current, chosen, frontier, accepted_edges_latest=None):
        """
        Chame isso dentro do loop. Grava 1 frame a cada sample_every.
        accepted_edges_latest: conjunto de arestas aceitas da folha mais recente, para destacar.
        """
        if not self.enabled or self.G is None:
            return
        self.step += 1
        if self.step % self.sample_every != 0:
            return
        if len(self.frames) >= self.max_frames:
            return

        if current is None:
            current_t = ()
        elif isinstance(current, (list, tuple, set)):
            current_t = tuple(int(x) for x in current)
        else:
            current_t = (int(current),)

        chosen_t = tuple(int(x) for x in chosen)
        frontier_t = tuple(int(x) for x in frontier)
        aes = None
        if accepted_edges_latest:
            aes = tuple(sorted(tuple(sorted(e)) for e in accepted_edges_latest))
        self.frames.append((current_t, chosen_t, frontier_t, aes))

    def end_anchor(self, name_prefix=None):
        """
        Renderiza e salva. Se tiver zero frames, não salva nada.
        """
        if not self.enabled or self.G is None:
            return None

        if not self.frames:
            return None

        name_prefix = name_prefix or f"anchor_{self.anchor}"
        out_path = self.out_dir / f"{name_prefix}.{self.fmt}"
        if self.fmt == "gif":
            if imageio is None:
                raise RuntimeError("imageio não disponível. Instale imageio para GIF ou use fmt='mp4'.")
            self._render_gif(out_path)
        elif self.fmt == "mp4":
            self._render_mp4(out_path)
        else:
            raise ValueError("fmt deve ser 'gif' ou 'mp4'.")

        saved = str(out_path)
        self.reset()
        return saved

    def _draw_frame(self, ax, current, chosen, frontier, accepted_edges_latest):
        # cores
        base_color = (0.6, 0.6, 0.6)
        chosen_color = (0.0, 0.6, 0.0)
        frontier_color = (0.9, 0.5, 0.0)
        current_color = (0.0, 0.4, 0.9)
        anchor_color = (0.9, 0.1, 0.1)

        # estado
        chosen_set = set(chosen)
        frontier_set = set(frontier)
        current_set = set(current)

        # desenha edges base
        nx.draw_networkx_edges(self.G, self.pos, ax=ax, width=1.0, alpha=self.edge_alpha)

        # nós em 4 grupos
        nodes_all = list(self.G.nodes())
        nodes_anchor = [self.anchor] if self.anchor in self.G else []
        nodes_current = [n for n in current_set if n not in nodes_anchor]
        nodes_chosen = [n for n in chosen_set if n not in nodes_anchor and n not in nodes_current]
        nodes_frontier = [n for n in frontier_set if n not in nodes_anchor and n not in nodes_current and n not in chosen_set]
        nodes_rest = [n for n in nodes_all if n not in chosen_set and n not in frontier_set and n not in nodes_current and n not in nodes_anchor]

        nx.draw_networkx_nodes(self.G, self.pos, nodelist=nodes_rest, node_size=40, node_color=[base_color], ax=ax)
        if nodes_chosen:
            nx.draw_networkx_nodes(self.G, self.pos, nodelist=nodes_chosen, node_size=60, node_color=[chosen_color], ax=ax)
        if nodes_frontier:
            nx.draw_networkx_nodes(self.G, self.pos, nodelist=nodes_frontier, node_size=70, node_color=[frontier_color], ax=ax)
        if nodes_current:
            nx.draw_networkx_nodes(self.G, self.pos, nodelist=nodes_current, node_size=90, node_color=[current_color], ax=ax)
        if nodes_anchor:
            nx.draw_networkx_nodes(self.G, self.pos, nodelist=nodes_anchor, node_size=120, node_color=[anchor_color], ax=ax)

        # label só do âncora
        if nodes_anchor:
            nx.draw_networkx_labels(self.G, self.pos, labels={self.anchor: str(self.anchor)}, font_size=9, ax=ax)

        # se tiver edges aceitas, destaca
        if accepted_edges_latest:
            # accepted_edges_latest é um iterável de pares (u, v)
            H = nx.Graph()
            for u, v in accepted_edges_latest:
                if self.G.has_node(u) and self.G.has_node(v) and self.G.has_edge(u, v):
                    H.add_edge(u, v)
            nx.draw_networkx_edges(H, self.pos, ax=ax, width=2.5, alpha=0.95)
        
        ax.set_title(f"anchor={self.anchor} current={list(current)} chosen={len(chosen)} frontier={len(frontier)}", fontsize=8)
        ax.set_axis_off()

    def _render_gif(self, out_path):
        imgs = []
        # renderiza em memória para GIF
        for current, chosen, frontier, accepted in self.frames:
            fig = plt.figure(figsize=(6, 6), dpi=self.dpi)
            ax = plt.gca()
            self._draw_frame(ax, current, chosen, frontier, accepted)
            fig.canvas.draw()
            img = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
            img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            imgs.append(img)
            plt.close(fig)
        imageio.mimsave(out_path, imgs, fps=self.fps)

    def _render_mp4(self, out_path):
        # usa matplotlib + FFMpegWriter
        from matplotlib.animation import FFMpegWriter
        metadata = dict(artist="TraversalTracer")
        writer = FFMpegWriter(fps=self.fps, metadata=metadata, bitrate=1800)
        fig = plt.figure(figsize=(6, 6), dpi=self.dpi)
        with writer.saving(fig, str(out_path), dpi=self.dpi):
            for current, chosen, frontier, accepted in self.frames:
                fig.clf()
                ax = fig.gca()
                self._draw_frame(ax, current, chosen, frontier, accepted)
                writer.grab_frame()
        plt.close(fig)

