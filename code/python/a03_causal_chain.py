# a03_causal_chain.py —— 因果链 / 结构相关（移植 R 10_analysis_03_causal_chain.R）
# H/α/载体 与心理-躯体结局的相关；phq9 对四 α 的多重回归（标准化系数）。
import numpy as np
import pandas as pd
import common as C

MODULE = "a03_causal_chain"
CORE = ["H_score", "alpha1", "alpha2", "alpha3", "alpha4",
        "carrier_a", "carrier_b", "carrier_c", "carrier_d"]
OUTC = ["phq9_total", "gad7_total", "ibs_symptom_count"]


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        # 相关阵
        cols = CORE + OUTC
        sub = d[cols].dropna()
        R = sub.corr().round(4)
        C.save_tab(R.reset_index().rename(columns={"index": "变量"}), "a03_相关阵", cohort, MODULE)
        # phq9 ~ α1-α4 标准化回归
        import statsmodels.api as sm
        reg = sub[["alpha1", "alpha2", "alpha3", "alpha4", "phq9_total"]].dropna()
        predictors = reg[["alpha1", "alpha2", "alpha3", "alpha4"]].astype(float)
        Xz = (predictors - predictors.mean()) / predictors.std(ddof=0)
        yz = (reg["phq9_total"] - reg["phq9_total"].mean()) / reg["phq9_total"].std(ddof=0)
        X = sm.add_constant(Xz)
        res = sm.OLS(yz, X).fit()
        coef = pd.DataFrame({
            "预测子": ["alpha1", "alpha2", "alpha3", "alpha4"],
            "standardized_beta": np.round(res.params.values[1:], 4),
            "p": np.round(res.pvalues.values[1:], 5),
        })
        coef.insert(0, "队列", cohort)
        C.save_tab(coef, "a03_phq9_对四α回归", cohort, MODULE)
        print("[a03] %s 完成 | phq9~α R²=%.3f" % (cohort, res.rsquared))
    print("[a03] 因果链/结构相关完成")


if __name__ == "__main__":
    main()
