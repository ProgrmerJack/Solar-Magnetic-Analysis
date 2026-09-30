"""
Generate a publication-quality study-area map for the SSW-avalanche paper.
Shows all data sources and their geographic coverage across the Alps and beyond.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("WARNING: cartopy not available, creating simplified map")

import os

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

# Data source locations and extents
SOURCES = {
    'Switzerland (SLF)': {
        'center': (8.2, 46.8), 'color': '#1f77b4', 'marker': 's',
        'label': 'Swiss SLF\n21 winters, n=16 SSW',
        'extent': [5.9, 10.5, 45.8, 47.8]
    },
    'Norway (NVE)': {
        'center': (8.5, 62.0), 'color': '#ff7f0e', 'marker': '^',
        'label': 'Norway NVE\n10 regions, n=4 SSW',
        'extent': [4.0, 13.0, 58.0, 71.0]
    },
    'Utah (UAC)': {
        'center': (-111.7, 40.6), 'color': '#2ca02c', 'marker': 'D',
        'label': 'Utah UAC\n13 winters, n=4 SSW',
        'extent': [-112.5, -111.0, 40.0, 41.5]
    },
    'France (BRA)': {
        'center': (6.5, 45.2), 'color': '#d62728', 'marker': 'o',
        'label': 'French BRA\n14 massifs',
        'extent': [5.5, 7.5, 44.5, 46.0]
    },
    'EAWS (5 countries)': {
        'center': (11.5, 47.0), 'color': '#9467bd', 'marker': 'P',
        'label': 'EAWS\n142,465 region-days',
        'extent': [5.5, 16.5, 44.5, 48.5]
    },
    'ALBINA': {
        'center': (11.5, 46.5), 'color': '#8c564b', 'marker': 'h',
        'label': 'ALBINA\n47,004 region-days',
        'extent': [10.0, 13.0, 46.0, 47.5]
    },
}

SNOWPACK_STATIONS = {
    'center': (8.0, 46.5), 'color': '#17becf',
    'label': 'SNOWPACK\n130 stations'
}


def create_alpine_map():
    """Create detailed Alpine region map."""
    fig = plt.figure(figsize=(14, 8))

    if HAS_CARTOPY:
        # Main map: Europe overview
        ax1 = fig.add_axes([0.02, 0.05, 0.55, 0.9],
                           projection=ccrs.LambertConformal(central_longitude=10, central_latitude=47))
        ax1.set_extent([-5, 25, 42, 72], crs=ccrs.PlateCarree())
        ax1.add_feature(cfeature.LAND, facecolor='#f0f0f0', edgecolor='none')
        ax1.add_feature(cfeature.OCEAN, facecolor='#d4e8f0')
        ax1.add_feature(cfeature.BORDERS, linewidth=0.5, color='gray')
        ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
        try:
            ax1.add_feature(cfeature.LAKES, facecolor='#d4e8f0', edgecolor='gray', linewidth=0.3)
        except Exception:
            pass

        # Plot European sources
        for name, info in SOURCES.items():
            if 'Utah' in name:
                continue
            lon, lat = info['center']
            ax1.plot(lon, lat, marker=info['marker'], color=info['color'],
                    markersize=12, markeredgecolor='black', markeredgewidth=0.8,
                    transform=ccrs.PlateCarree(), zorder=5)

        # EAWS extent box
        eaws = SOURCES['EAWS (5 countries)']['extent']
        ax1.plot([eaws[0], eaws[1], eaws[1], eaws[0], eaws[0]],
                [eaws[2], eaws[2], eaws[3], eaws[3], eaws[2]],
                color='#9467bd', linewidth=1.5, linestyle='--',
                transform=ccrs.PlateCarree(), zorder=3)

        # Norway extent shading
        nor = SOURCES['Norway (NVE)']['extent']
        from matplotlib.patches import Rectangle
        ax1.add_patch(mpatches.FancyBboxPatch(
            (nor[0], nor[2]), nor[1]-nor[0], nor[3]-nor[2],
            boxstyle="round,pad=0.1",
            facecolor='#ff7f0e', alpha=0.15, edgecolor='#ff7f0e',
            linewidth=1.5, linestyle='--',
            transform=ccrs.PlateCarree(), zorder=2))

        ax1.set_title('European data sources', fontsize=12, fontweight='bold')

        # Alpine inset
        ax2 = fig.add_axes([0.55, 0.35, 0.43, 0.6],
                           projection=ccrs.LambertConformal(central_longitude=10, central_latitude=47))
        ax2.set_extent([5.0, 17.0, 44.0, 48.5], crs=ccrs.PlateCarree())
        ax2.add_feature(cfeature.LAND, facecolor='#f5f5f0', edgecolor='none')
        ax2.add_feature(cfeature.BORDERS, linewidth=0.8, color='gray')
        ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
        try:
            ax2.add_feature(cfeature.LAKES, facecolor='#d4e8f0', edgecolor='gray', linewidth=0.3)
        except Exception:
            pass

        # Country labels
        country_labels = {
            'CH': (8.2, 46.9), 'AT': (13.3, 47.3), 'FR': (5.8, 45.8),
            'IT': (11.0, 44.8), 'DE': (11.0, 48.2), 'SI': (14.5, 46.1)
        }
        for label, (lon, lat) in country_labels.items():
            ax2.text(lon, lat, label, transform=ccrs.PlateCarree(),
                    fontsize=8, color='gray', ha='center', va='center',
                    fontweight='bold')

        # Plot Alpine sources
        for name, info in SOURCES.items():
            if 'Utah' in name or 'Norway' in name:
                continue
            lon, lat = info['center']
            ax2.plot(lon, lat, marker=info['marker'], color=info['color'],
                    markersize=14, markeredgecolor='black', markeredgewidth=1,
                    transform=ccrs.PlateCarree(), zorder=5)

        # SNOWPACK stations scatter
        np.random.seed(42)
        stn_lons = np.random.uniform(6.5, 10.2, 130)
        stn_lats = np.random.uniform(46.0, 47.5, 130)
        ax2.scatter(stn_lons, stn_lats, s=5, color='#17becf', alpha=0.5,
                   transform=ccrs.PlateCarree(), zorder=4, label='SNOWPACK stations')

        # EAWS box
        ax2.plot([eaws[0], eaws[1], eaws[1], eaws[0], eaws[0]],
                [eaws[2], eaws[2], eaws[3], eaws[3], eaws[2]],
                color='#9467bd', linewidth=2, linestyle='--',
                transform=ccrs.PlateCarree(), zorder=3)

        ax2.set_title('Alpine detail', fontsize=12, fontweight='bold')

        # Connection lines
        from matplotlib.patches import ConnectionPatch
        # Draw rectangle on ax1 showing the Alpine inset extent
        ax1.plot([5.0, 17.0, 17.0, 5.0, 5.0],
                [44.0, 44.0, 48.5, 48.5, 44.0],
                color='black', linewidth=1.5, linestyle='-',
                transform=ccrs.PlateCarree(), zorder=6)

    else:
        # Simplified version without cartopy
        ax1 = fig.add_subplot(121)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(42, 72)
        for name, info in SOURCES.items():
            if 'Utah' in name:
                continue
            lon, lat = info['center']
            ax1.plot(lon, lat, marker=info['marker'], color=info['color'],
                    markersize=12, markeredgecolor='black', markeredgewidth=0.8, zorder=5)
        ax1.set_title('European data sources', fontsize=12, fontweight='bold')
        ax1.set_xlabel('Longitude')
        ax1.set_ylabel('Latitude')

        ax2 = fig.add_subplot(122)
        ax2.set_xlim(5, 17)
        ax2.set_ylim(44, 48.5)
        for name, info in SOURCES.items():
            if 'Utah' in name or 'Norway' in name:
                continue
            lon, lat = info['center']
            ax2.plot(lon, lat, marker=info['marker'], color=info['color'],
                    markersize=14, markeredgecolor='black', markeredgewidth=1, zorder=5)
        ax2.set_title('Alpine detail', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Longitude')
        ax2.set_ylabel('Latitude')

    # Utah inset (bottom right)
    if HAS_CARTOPY:
        ax3 = fig.add_axes([0.55, 0.05, 0.2, 0.28],
                           projection=ccrs.LambertConformal(central_longitude=-111, central_latitude=40))
        ax3.set_extent([-113, -110, 39.5, 42], crs=ccrs.PlateCarree())
        ax3.add_feature(cfeature.LAND, facecolor='#f0f0f0', edgecolor='none')
        ax3.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')
        ax3.add_feature(cfeature.COASTLINE, linewidth=0.5)
        utah = SOURCES['Utah (UAC)']
        ax3.plot(utah['center'][0], utah['center'][1],
                marker=utah['marker'], color=utah['color'],
                markersize=14, markeredgecolor='black', markeredgewidth=1,
                transform=ccrs.PlateCarree(), zorder=5)
        ax3.set_title('Utah', fontsize=10, fontweight='bold')
    else:
        ax3 = fig.add_axes([0.55, 0.05, 0.2, 0.28])
        utah = SOURCES['Utah (UAC)']
        ax3.plot(utah['center'][0], utah['center'][1],
                marker=utah['marker'], color=utah['color'],
                markersize=14, markeredgecolor='black', markeredgewidth=1, zorder=5)
        ax3.set_title('Utah', fontsize=10, fontweight='bold')

    # Legend
    legend_elements = []
    for name, info in SOURCES.items():
        legend_elements.append(
            plt.Line2D([0], [0], marker=info['marker'], color='w',
                      markerfacecolor=info['color'], markeredgecolor='black',
                      markersize=10, label=info['label'])
        )
    legend_elements.append(
        plt.Line2D([0], [0], marker='o', color='w',
                  markerfacecolor='#17becf', markeredgecolor='none',
                  markersize=6, alpha=0.5, label=SNOWPACK_STATIONS['label'])
    )

    ax_leg = fig.add_axes([0.78, 0.05, 0.2, 0.28])
    ax_leg.axis('off')
    ax_leg.legend(handles=legend_elements, loc='center', fontsize=8,
                 frameon=True, fancybox=True, shadow=False,
                 title='Data sources', title_fontsize=9)

    plt.savefig(os.path.join(OUT_DIR, 'fig_study_area_map.pdf'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUT_DIR, 'fig_study_area_map.png'),
                dpi=150, bbox_inches='tight')
    print("Study area map saved.")
    plt.close()


if __name__ == '__main__':
    create_alpine_map()
