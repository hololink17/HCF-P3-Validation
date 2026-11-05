# a16_hierarchical_regression.py —— 层级回归（移植 R 10_analysis_16_hierarchical_regression.R）
# 预测抑郁(PHQ≥10)：M0=人口学；M1=+α1-α4；M2=+四载体；M3=+H_score。
# 报告 R² / ΔR² / F / p（两队列）。
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import common as C

MODULE = "a16_hierarchical_regression"
DEMO = ["age_group", "gender", "education"]


def dummies(d):
    X = pd.DataFrame(index=d.index)
    for c in DEMO:
        if c in d.columns:
            X = pd.concat([X, pd.get_dummies(d[c].astype(str), prefix=c, drop_first=True)],
                          axis=1)
    return X


def run(d):
    out = []
    required = (["depression"] + DEMO + ["alpha1", "alpha2", "alpha3", "alpha4"] +
                ["carrier_a", "carrier_b", "carrier_c", "carrier_d", "H_score"])
    # 所有嵌套模型固定同一完整病例样本，避免 ΔR² 混入样本构成变化。
    base = d[required].dropna()
    y = base["depression"]
    X0 = dummies(base)
    X1 = pd.concat([X0, base[["alpha1", "alpha2", "alpha3", "alpha4"]].astype(float)], axis=1)
    X2 = pd.concat([X1, base[["carrier_a", "carrier_b", "carrier_c", "carrier_d"]].astype(float)], axis=1)
    X3 = pd.concat([X2, base[["H_score"]].astype(float)], axis=1)
    models = [("M0_人口学", X0), ("M1_+α", X1), ("M2_+载体", X2), ("M3_+H", X3)]
    prev_r2 = 0.0
    for name, X in models:
        Xc = sm.add_constant(X.dropna(axis=1, how="all"))
        dd = Xc.dropna()
        yy = y.loc[dd.index].astype(float)
        if yy.nunique() < 2 or dd.shape[1] < 1:
            continue
        # 同一二项最大似然模型使用 GLM/IRLS，避免 Newton 数值容差造成假未收敛。
        res = sm.GLM(yy, dd, family=sm.families.Binomial()).fit()
        if not bool(res.converged):
            raise RuntimeError("%s Logistic模型未收敛" % name)
        r2 = 1-res.llf/res.llnull
        overall_p = stats.chi2.sf(2*(res.llf-res.llnull), res.df_model)
        out.append((name, round(overall_p, 5), round(r2, 4),
                    round(r2 - prev_r2, 4), dd.shape[1], int(len(yy))))
        prev_r2 = r2
    return pd.DataFrame(out, columns=["模型", "整体p", "McFadden_R²", "ΔR²", "预测子数", "n"])


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        C.save_tab(run(d), "a16_层级回归_抑郁", cohort, MODULE)
        print("[a16] %s 层级回归完成" % cohort)
    print("[a16] 层级回归完成")


if __name__ == "__main__":
    main()
