"""DOY-matched mediation analysis — fixes seasonal confounding."""
import pandas as pd, numpy as np
from scipy import stats
from numpy.linalg import lstsq
np.random.seed(42)

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_cat.index = ssw_cat.index.tz_localize(None)
ssw_dates = ssw_cat[ssw_cat['type'] == 'M'].index
ssw_dates = ssw_dates[(ssw_dates >= '1998-01-01') & (ssw_dates <= '2019-12-31')]

disp_dates = pd.to_datetime(['1998-12-15', '1999-02-26', '2001-02-11', '2004-01-05', '2006-01-21',
                              '2007-02-24', '2008-02-22', '2010-02-09', '2012-01-11', '2019-01-01'])

merged = panel.join(era5[['t2m_K', 'sf_mm', 'wind_speed']], how='inner')
merged = merged[merged['is_winter'] == 1]

# Build DOY-matched treatment/control pairs
treat_rows = []
ctrl_rows = []

for ssw_date in ssw_dates:
    ssw_window = pd.date_range(ssw_date, ssw_date + pd.Timedelta(days=15))
    ssw_type = 'displacement' if ssw_date in disp_dates else 'split'
    
    for day in ssw_window:
        if day in merged.index:
            row = merged.loc[day]
            treat_rows.append({
                'date': day, 'doy': day.dayofyear,
                't2m_K': row['t2m_K'], 'sf_mm': row.get('sf_mm', np.nan),
                'wind_speed': row.get('wind_speed', np.nan),
                'aval': row['aai_all_dry'], 'ssw_date': ssw_date, 'type': ssw_type,
            })
    
    # Control: same DOYs from non-SSW years
    for year in range(1998, 2020):
        if year == ssw_date.year:
            continue
        for doy in [d.dayofyear for d in ssw_window]:
            try:
                ctrl_date = pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=doy-1)
            except:
                continue
            if ctrl_date not in merged.index:
                continue
            near_ssw = any(abs((ctrl_date - s).days) < 30 for s in ssw_dates)
            if near_ssw:
                continue
            row = merged.loc[ctrl_date]
            ctrl_rows.append({
                'date': ctrl_date, 'doy': ctrl_date.dayofyear,
                't2m_K': row['t2m_K'], 'sf_mm': row.get('sf_mm', np.nan),
                'wind_speed': row.get('wind_speed', np.nan),
                'aval': row['aai_all_dry'], 'ssw_date': ssw_date, 'type': ssw_type,
            })

treat_df = pd.DataFrame(treat_rows)
ctrl_df = pd.DataFrame(ctrl_rows)
print("Treatment days: %d, Control days: %d" % (len(treat_df), len(ctrl_df)))

# DOY-matched comparison
print("\n=== DOY-MATCHED SSW vs CONTROL ===")
print("Treatment aval rate: %.3f" % treat_df['aval'].mean())
print("Control aval rate: %.3f" % ctrl_df['aval'].mean())
rr = treat_df['aval'].mean() / ctrl_df['aval'].mean()
print("Rate ratio: %.3f" % rr)
u, p = stats.mannwhitneyu(ctrl_df['aval'], treat_df['aval'], alternative='two-sided')
print("Mann-Whitney P = %.6f" % p)

dt = treat_df['t2m_K'].mean() - ctrl_df['t2m_K'].mean()
print("\nDelta T: %+.3f K" % dt)
t, p_t = stats.ttest_ind(treat_df['t2m_K'], ctrl_df['t2m_K'])
print("T-test for T: P = %.6f" % p_t)

dsf = treat_df['sf_mm'].mean() - ctrl_df['sf_mm'].mean()
print("Delta snowfall: %+.4f mm/d" % dsf)

dw = treat_df['wind_speed'].mean() - ctrl_df['wind_speed'].mean()
print("Delta wind: %+.4f m/s" % dw)

# Mediation analysis
print("\n=== DOY-MATCHED MEDIATION (Baron & Kenny) ===")
treat_df['treatment'] = 1
ctrl_df['treatment'] = 0
combined = pd.concat([treat_df, ctrl_df], ignore_index=True).dropna()

# Step 1: X -> Y
X1 = np.column_stack([combined['treatment'].values, np.ones(len(combined))])
y1 = combined['aval'].values
b1, _, _, _ = lstsq(X1, y1, rcond=None)
resid1 = y1 - X1 @ b1
mse1 = np.sum(resid1**2) / (len(y1) - 2)
se1 = np.sqrt(mse1 * np.linalg.inv(X1.T @ X1).diagonal())
t_total = b1[0] / se1[0]
p_total = 2 * (1 - stats.t.cdf(abs(t_total), len(y1)-2))
print("Step 1 (X->Y): b=%.4f, SE=%.4f, t=%.3f, P=%.6f" % (b1[0], se1[0], t_total, p_total))

# Step 2: X -> M (temperature)
m2 = combined['t2m_K'].values
b2, _, _, _ = lstsq(X1, m2, rcond=None)
resid2 = m2 - X1 @ b2
mse2 = np.sum(resid2**2) / (len(m2) - 2)
se2 = np.sqrt(mse2 * np.linalg.inv(X1.T @ X1).diagonal())
t_xm = b2[0] / se2[0]
p_xm = 2 * (1 - stats.t.cdf(abs(t_xm), len(m2)-2))
print("Step 2 (X->M): b=%.4f, SE=%.4f, t=%.3f, P=%.6f" % (b2[0], se2[0], t_xm, p_xm))

# Step 3: X + M -> Y
X3 = np.column_stack([combined['treatment'].values, combined['t2m_K'].values, np.ones(len(combined))])
b3, _, _, _ = lstsq(X3, y1, rcond=None)
resid3 = y1 - X3 @ b3
mse3 = np.sum(resid3**2) / (len(y1) - 3)
se3 = np.sqrt(mse3 * np.linalg.inv(X3.T @ X3).diagonal())

a = b2[0]  # X -> M
b = b3[1]  # M -> Y | X
indirect = a * b
direct = b3[0]

sobel_se = np.sqrt(a**2 * se3[1]**2 + b**2 * se2[0]**2)
sobel_z = indirect / sobel_se if sobel_se > 0 else 0
sobel_p = 2 * (1 - stats.norm.cdf(abs(sobel_z)))
prop = indirect / b1[0] if abs(b1[0]) > 0.001 else float('nan')

print("Direct effect (c'): %.4f" % direct)
print("Mediator effect (b): %.4f" % b)
print("Indirect (a*b): %.4f" % indirect)
print("Sobel z=%.3f, P=%.6f" % (sobel_z, sobel_p))
print("Proportion mediated: %.1f%%" % (prop * 100))

# Bootstrap CI
n_boot = 5000
boot_indirect = []
for _ in range(n_boot):
    idx = np.random.choice(len(combined), len(combined), replace=True)
    bdf = combined.iloc[idx]
    bX = np.column_stack([bdf['treatment'].values, np.ones(len(bdf))])
    bm = bdf['t2m_K'].values
    bb2, _, _, _ = lstsq(bX, bm, rcond=None)
    bX3 = np.column_stack([bdf['treatment'].values, bdf['t2m_K'].values, np.ones(len(bdf))])
    by = bdf['aval'].values
    try:
        bb3, _, _, _ = lstsq(bX3, by, rcond=None)
        boot_indirect.append(bb2[0] * bb3[1])
    except:
        pass

ci = np.percentile(boot_indirect, [2.5, 97.5])
print("Bootstrap 95%% CI: [%.4f, %.4f]" % (ci[0], ci[1]))

# Wind mediation
print("\n=== WIND MEDIATION ===")
m2w = combined['wind_speed'].values
b2w, _, _, _ = lstsq(X1, m2w, rcond=None)
resid2w = m2w - X1 @ b2w
mse2w = np.sum(resid2w**2) / (len(m2w) - 2)
se2w = np.sqrt(mse2w * np.linalg.inv(X1.T @ X1).diagonal())
print("Path a (X->Wind): b=%.4f, P=%.6f" % (b2w[0], 2*(1-stats.t.cdf(abs(b2w[0]/se2w[0]), len(m2w)-2))))

X3w = np.column_stack([combined['treatment'].values, combined['wind_speed'].values, np.ones(len(combined))])
b3w, _, _, _ = lstsq(X3w, y1, rcond=None)
indirect_w = b2w[0] * b3w[1]
print("Path b (Wind->Y|X): b=%.4f" % b3w[1])
print("Indirect (wind): %.4f" % indirect_w)
print("Proportion: %.1f%%" % (indirect_w / b1[0] * 100 if abs(b1[0]) > 0.001 else 0))

# Snowfall mediation
print("\n=== SNOWFALL MEDIATION ===")
m2s = combined['sf_mm'].values
b2s, _, _, _ = lstsq(X1, m2s, rcond=None)
resid2s = m2s - X1 @ b2s
mse2s = np.sum(resid2s**2) / (len(m2s) - 2)
se2s = np.sqrt(mse2s * np.linalg.inv(X1.T @ X1).diagonal())
print("Path a (X->SF): b=%.4f, P=%.6f" % (b2s[0], 2*(1-stats.t.cdf(abs(b2s[0]/se2s[0]), len(m2s)-2))))

X3s = np.column_stack([combined['treatment'].values, combined['sf_mm'].values, np.ones(len(combined))])
b3s, _, _, _ = lstsq(X3s, y1, rcond=None)
indirect_s = b2s[0] * b3s[1]
print("Path b (SF->Y|X): b=%.4f" % b3s[1])
print("Indirect (snowfall): %.4f" % indirect_s)
print("Proportion: %.1f%%" % (indirect_s / b1[0] * 100 if abs(b1[0]) > 0.001 else 0))

# DISPLACEMENT ONLY
print("\n" + "="*60)
print("=== DISPLACEMENT-ONLY DOY-MATCHED ===")
print("="*60)
td = treat_df[treat_df['type'] == 'displacement']
cd = ctrl_df[ctrl_df['type'] == 'displacement']
print("Treatment: %d, Control: %d" % (len(td), len(cd)))
print("Treat aval: %.3f, Ctrl aval: %.3f" % (td['aval'].mean(), cd['aval'].mean()))
rr_d = td['aval'].mean() / cd['aval'].mean()
print("RR = %.3f" % rr_d)
u_d, p_d = stats.mannwhitneyu(cd['aval'], td['aval'], alternative='two-sided')
print("MW P = %.6f" % p_d)
dt_d = td['t2m_K'].mean() - cd['t2m_K'].mean()
print("Delta T: %+.3f K" % dt_d)
t_d, p_td = stats.ttest_ind(td['t2m_K'], cd['t2m_K'])
print("T P = %.6f" % p_td)

# Displacement mediation
td2 = td.copy()
cd2 = cd.copy()
td2['treatment'] = 1
cd2['treatment'] = 0
comb_d = pd.concat([td2, cd2], ignore_index=True).dropna()

X1d = np.column_stack([comb_d['treatment'].values, np.ones(len(comb_d))])
y1d = comb_d['aval'].values
b1d, _, _, _ = lstsq(X1d, y1d, rcond=None)

m2d = comb_d['t2m_K'].values
b2d, _, _, _ = lstsq(X1d, m2d, rcond=None)
resid2d = m2d - X1d @ b2d
mse2d = np.sum(resid2d**2) / (len(m2d) - 2)
se2d = np.sqrt(mse2d * np.linalg.inv(X1d.T @ X1d).diagonal())

X3d = np.column_stack([comb_d['treatment'].values, comb_d['t2m_K'].values, np.ones(len(comb_d))])
b3d, _, _, _ = lstsq(X3d, y1d, rcond=None)
resid3d = y1d - X3d @ b3d
mse3d = np.sum(resid3d**2) / (len(y1d) - 3)
se3d = np.sqrt(mse3d * np.linalg.inv(X3d.T @ X3d).diagonal())

a_d = b2d[0]
b_d = b3d[1]
indirect_d = a_d * b_d

sobel_se_d = np.sqrt(a_d**2 * se3d[1]**2 + b_d**2 * se2d[0]**2)
sobel_z_d = indirect_d / sobel_se_d if sobel_se_d > 0 else 0
sobel_p_d = 2 * (1 - stats.norm.cdf(abs(sobel_z_d)))

print("\nMediation (displacement):")
print("Path a (X->T): %.4f, P=%.6f" % (a_d, 2*(1-stats.t.cdf(abs(a_d/se2d[0]), len(m2d)-2))))
print("Path b (T->Y|X): %.4f" % b_d)
print("Indirect: %.4f" % indirect_d)
print("Sobel z=%.3f, P=%.6f" % (sobel_z_d, sobel_p_d))
if abs(b1d[0]) > 0.001:
    print("Proportion mediated: %.1f%%" % (indirect_d / b1d[0] * 100))

# Bootstrap
boot_ind_d = []
for _ in range(5000):
    idx = np.random.choice(len(comb_d), len(comb_d), replace=True)
    bdf = comb_d.iloc[idx]
    bX = np.column_stack([bdf['treatment'].values, np.ones(len(bdf))])
    bm = bdf['t2m_K'].values
    bb2, _, _, _ = lstsq(bX, bm, rcond=None)
    bX3 = np.column_stack([bdf['treatment'].values, bdf['t2m_K'].values, np.ones(len(bdf))])
    by = bdf['aval'].values
    try:
        bb3, _, _, _ = lstsq(bX3, by, rcond=None)
        boot_ind_d.append(bb2[0] * bb3[1])
    except:
        pass
ci_d = np.percentile(boot_ind_d, [2.5, 97.5])
print("Bootstrap 95%% CI: [%.4f, %.4f]" % (ci_d[0], ci_d[1]))

# Multi-mediator model (T + wind + snowfall simultaneously)
print("\n=== MULTI-MEDIATOR (T + Wind + SF) ===")
X_multi = np.column_stack([
    combined['treatment'].values,
    combined['t2m_K'].values,
    combined['wind_speed'].values,
    combined['sf_mm'].values,
    np.ones(len(combined))
])
b_multi, _, _, _ = lstsq(X_multi, y1, rcond=None)
resid_m = y1 - X_multi @ b_multi
ss_res = np.sum(resid_m**2)
ss_tot = np.sum((y1 - y1.mean())**2)
r2_multi = 1 - ss_res / ss_tot

print("Direct (treatment|T,Wind,SF): %.4f" % b_multi[0])
print("Temperature: %.4f" % b_multi[1])
print("Wind: %.4f" % b_multi[2])
print("Snowfall: %.4f" % b_multi[3])
print("R-squared: %.4f" % r2_multi)

total_indirect = b1[0] - b_multi[0]
print("Total indirect (all mediators): %.4f" % total_indirect)
print("Total proportion mediated: %.1f%%" % (total_indirect / b1[0] * 100 if abs(b1[0]) > 0.001 else 0))

print("\nDONE")
