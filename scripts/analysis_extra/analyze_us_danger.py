import pandas as pd
import numpy as np
from scipy import stats

# Analyze US danger ratings during 2021 SSW vs 2022 control
ssw_df = pd.read_csv(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\us_danger_ratings_2021_ssw.csv')
ctrl_df = pd.read_csv(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\us_danger_ratings_2022_control.csv')

# Clean up dates, drop NaN danger ratings
ssw_df['date'] = pd.to_datetime(ssw_df['date'])
ctrl_df['date'] = pd.to_datetime(ctrl_df['date'])
ssw_df = ssw_df.dropna(subset=['danger_rating'])
ctrl_df = ctrl_df.dropna(subset=['danger_rating'])

# SSW onset: 2021-01-05
# Define windows: pre-SSW (-30 to -1), post-SSW (0 to +30)
ssw_onset = pd.Timestamp('2021-01-05')
ssw_df['days_from_onset'] = (ssw_df['date'] - ssw_onset).dt.days

# Narrow to the relevant period
ssw_period = ssw_df[(ssw_df['days_from_onset'] >= -30) & (ssw_df['days_from_onset'] <= 30)]

# Daily mean danger across all centers
daily = ssw_period.groupby('days_from_onset')['danger_rating'].agg(['mean','std','count'])
print('=== Daily mean danger rating around 2021-01-05 SSW ===')
for d in range(-15, 16):
    if d in daily.index:
        row = daily.loc[d]
        print(f'Day {d:+3d}: mean={row["mean"]:.2f}, n={int(row["count"])}')

# Pre vs Post SSW comparison
pre = ssw_period[ssw_period['days_from_onset'].between(-15, -1)]['danger_rating']
post = ssw_period[ssw_period['days_from_onset'].between(0, 15)]['danger_rating']
print(f'\nPre-SSW mean danger: {pre.mean():.3f} +/- {pre.std():.3f} (n={len(pre)})')
print(f'Post-SSW mean danger: {post.mean():.3f} +/- {post.std():.3f} (n={len(post)})')
t_stat, t_p = stats.ttest_ind(pre, post)
u_stat, u_p = stats.mannwhitneyu(pre, post, alternative='two-sided')
print(f'T-test: t={t_stat:.3f}, P={t_p:.4f}')
print(f'Mann-Whitney: P={u_p:.4f}')

# Same comparison for 2022 control (no SSW)
ctrl_onset = pd.Timestamp('2022-01-05')
ctrl_df['days_from_onset'] = (ctrl_df['date'] - ctrl_onset).dt.days
ctrl_period = ctrl_df[(ctrl_df['days_from_onset'] >= -30) & (ctrl_df['days_from_onset'] <= 30)]
ctrl_pre = ctrl_period[ctrl_period['days_from_onset'].between(-15, -1)]['danger_rating']
ctrl_post = ctrl_period[ctrl_period['days_from_onset'].between(0, 15)]['danger_rating']
print(f'\n2022 Control:')
print(f'Pre mean: {ctrl_pre.mean():.3f}, Post mean: {ctrl_post.mean():.3f}')
t2, p2 = stats.ttest_ind(ctrl_pre, ctrl_post)
print(f'T-test: t={t2:.3f}, P={p2:.4f}')

# Per-center analysis
print('\n=== Per-center danger change (SSW period) ===')
centers = ssw_period['center'].unique()
decreases = 0
for c in sorted(centers):
    c_data = ssw_period[ssw_period['center'] == c]
    c_pre = c_data[c_data['days_from_onset'].between(-15, -1)]['danger_rating'].mean()
    c_post = c_data[c_data['days_from_onset'].between(0, 15)]['danger_rating'].mean()
    change = c_post - c_pre
    marker = 'v' if change < 0 else '^'
    if change < 0:
        decreases += 1
    print(f'{c[:30]:30s}: pre={c_pre:.2f} post={c_post:.2f} change={change:+.2f} {marker}')

print(f'\nDecreases: {decreases}/{len(centers)}')
sign_p = stats.binom_test(decreases, len(centers), 0.5) if hasattr(stats, 'binom_test') else 'N/A'
print(f'Sign test P={sign_p}')
