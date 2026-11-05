# =============================================================
# 00_shared_rules.R —— 跨平台共用规则层（唯一真值）
# 供 01_clean_*.R / 02_scoring.R / 10_analysis_*.R 调用；两个平台只认这里的口径
# 依据：docs/scoring_rules_v1_20260927.md + docs/分析用数据字典_v1_20260927.md
# 2026-09-27 用户裁定：
#   ① 疾病题 "不知道" → 记 0（同 Prolific：派生阳性数 = 原始值 1 的计数，2 全部落入阴性，已实证）
#   ② 年龄统一分箱为 18-35 / 36-45 / 46-60 / 60+
#   ③ 见数 A/B 两组不分层（同平台、同问卷、投放时间随机）
#   ④ 见数睡眠障碍变量作废 → 改用 sleep_symptom_count（0–3）
# 2026-09-28 增补：
#   ⑤ 版本标签统一为这里的两个常量（两队列一致，见下）
# =============================================================

# =============================================================
# 版本标签（两队列一致；仅作版本标识，不含内部规则）
# =============================================================
Q_VERSION    <- "科研-第一套题"
ALGO_VERSION <- "算法v1.1"

VERSION_TAG <- sprintf("问卷与映射 %s ｜ 算法参数 %s", Q_VERSION, ALGO_VERSION)

# ---------- 通用工具 ----------
lead_num <- function(x) suppressWarnings(as.numeric(str_extract(str_trim(as.character(x)), "^\\d+")))
has_code <- function(x, code) str_detect(as.character(x), paste0("(?<![0-9])", code, "(?![0-9])"))

# ---------- 年龄：统一分箱 ----------
AGE_BREAKS <- c(-Inf, 35, 45, 60, Inf)
AGE_LABELS <- c("18-35", "36-45", "46-60", "60+")

age_to_group <- function(age_num) {
  # 数值年龄（Prolific）→ 统一四档
  cut(as.numeric(age_num), breaks = AGE_BREAKS, labels = AGE_LABELS, right = TRUE)
}
age_label_to_group <- function(x) {
  # 见数已是四档文本，仅做规范化（防止全角/空格差异）
  x <- str_trim(as.character(x))
  x <- str_replace_all(x, "—|–|－", "-")
  factor(x, levels = AGE_LABELS)
}

# ---------- 性别对照（Prolific 1/2/3 ↔ 见数 男/女/不愿说）----------
# 实证：gender 1 → Sex=Male 251/254；2 → Female 206/212；3 → Prefer not to say（n=4）
GENDER_MAP <- c("1" = "男", "2" = "女", "3" = "不愿说")
norm_gender <- function(x) {
  x <- str_trim(as.character(x))
  out <- GENDER_MAP[x]
  out[is.na(out)] <- x[is.na(out)]
  factor(out, levels = c("男", "女", "不愿说"))
}

# ---------- 教育对照（4 档，两边选项一一对应）----------
# 见数：初中及以下 / 高中·中专 / 大专·本科 / 硕士及以上
# Prolific：1 Middle school or below / 2 High school·Vocational·Associate /
#           3 Bachelor's / 4 Master's or above
EDU_LEVELS <- c("初中及以下", "高中/中专", "大专/本科", "硕士及以上")
EDU_MAP_PROLIFIC <- c("1" = "初中及以下", "2" = "高中/中专", "3" = "大专/本科", "4" = "硕士及以上")
norm_education <- function(x, from = c("jianshu", "prolific")) {
  from <- match.arg(from)
  x <- str_trim(as.character(x))
  out <- if (from == "prolific") EDU_MAP_PROLIFIC[x] else
    dplyr::recode(x, "高中/中专" = "高中/中专", "大专/本科" = "大专/本科",
                  "初中及以下" = "初中及以下", "硕士及以上" = "硕士及以上")
  factor(out, levels = EDU_LEVELS)
}

# ---------- 怀孕（Prolific 0/1/2 ↔ 见数 否/是/不知道）----------
PREG_MAP <- c("0" = "否", "1" = "是", "2" = "不知道")
norm_pregnant <- function(x, from = c("jianshu", "prolific")) {
  from <- match.arg(from)
  x <- str_trim(as.character(x))
  out <- if (from == "prolific") PREG_MAP[x] else x
  factor(out, levels = c("是", "否", "不知道"))
}

# ---------- 疾病题：是=1，否/不知道/以上都没有=0 ----------
# 用户裁定：参考 Prolific（"不知道"记 0）
disease_from_text <- function(x) {
  x <- str_trim(as.character(x))
  dplyr::case_when(
    str_detect(x, "^1\\.|^是$|^1$") ~ 1,
    is.na(x) | x == "" ~ NA_real_,
    TRUE ~ 0                      # 2.否 / 3.不知道 / 4.以上都没有 全部记 0
  )
}

# ---------- IBS：选"腹痛或腹胀" 且 至少一项(腹泻/便秘/交替) ----------
ibs_positive_from_code <- function(x) {
  ifelse(is.na(x), NA_real_,
         as.numeric(has_code(x, "1") & (has_code(x, "2") | has_code(x, "3") | has_code(x, "4"))))
}

# ---------- 睡眠：症状计数 0–3（见数无 0 值，仅作严重度）----------
# 见数 A 组：单选文本；B 组：多选子列（"1|2"）
# Prolific：数字编码 "1;2;3"，4 = 以上都没有
sleep_count_from_text <- function(x) {
  x <- str_trim(as.character(x))
  ifelse(is.na(x) | x == "", NA_real_,
         as.numeric(str_detect(x, "入睡困难")) + as.numeric(str_detect(x, "夜间频繁醒来")) +
           as.numeric(str_detect(x, "早醒")))
}
sleep_count_from_code <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x) | x == "", NA_real_,
         as.numeric(has_code(x, "1")) + as.numeric(has_code(x, "2")) + as.numeric(has_code(x, "3")))
}

# ---------- PHQ-9 / GAD-7：前缀数字 − 1 ----------
likert0_3 <- function(x) { v <- lead_num(x); ifelse(is.na(v), NA_real_, v - 1) }

# 公开分析直接读取数值锚点；本文件不包含作答到分数的转换规则。

# =============================================================
# 分型规则（唯一真值）—— 2026-09-28 复核
# 依据：Prolific 原始 代码1「三波数据_数据分析-数据清洗合并」逐字定义
#   ① 每个 α 各自按本队列 P33 / P67 分「低 / 中 / 高」
#        x <= P33 → 低 ；x >= P67 → 高 ；否则 → 中
#   ② profile = α1cat-α2cat-α3cat-α4cat
#   ③ 四者同档 → 归类，否则「其他」
#        低-低-低-低 → 耗竭型 ｜ 中-中-中-中 → 均衡型 ｜ 高-高-高-高 → 高意识型
#   ④ self_ok     = phq9_total < 10
#      group_type = 耗竭型&正常→潜在抑郁 / 耗竭型&异常→典型抑郁 / 高意识型&正常→真正健康 / 其余 NA
#      possible_depression（潜在痛苦） = self_ok & alpha_group == "耗竭型"
# 交叉验证：本实现复现 Prolific 参考 alpha_group 471/471 逐字一致，
#           分组 H 均值 0.106 / 0.203 / 0.363 与参考完全吻合（2026-09-28）
# ⚠ 禁止改用「H 值三分位」或「α₁ 三分位」——那会得出 154/154/154 且无「其他」，
#    与参考口径不符（详见 docs/scoring_rules 与本次审计记录）
# =============================================================
ALPHA_AC <- c("alpha1", "alpha2", "alpha3", "alpha4")

add_alpha_typing <- function(d, ac = ALPHA_AC) {
  stopifnot(all(ac %in% names(d)))
  for (a in ac) {
    q <- stats::quantile(d[[a]], c(1 / 3, 2 / 3), na.rm = TRUE)
    d[[paste0(a, "_cat")]] <- ifelse(d[[a]] <= q[1], "低",
                              ifelse(d[[a]] >= q[2], "高", "中"))
  }
  d$profile <- paste(d$alpha1_cat, d$alpha2_cat, d$alpha3_cat, d$alpha4_cat, sep = "-")
  complete_alpha <- stats::complete.cases(d[, paste0(ac, "_cat"), drop = FALSE])
  d$profile[!complete_alpha] <- NA_character_
  d$alpha_group <- dplyr::case_when(
    !complete_alpha ~ NA_character_,
    d$profile == "低-低-低-低" ~ "耗竭型",
    d$profile == "中-中-中-中" ~ "均衡型",
    d$profile == "高-高-高-高" ~ "高意识型",
    TRUE ~ "其他"
  )
  d$self_ok        <- d$phq9_total < 10
  d$group_type     <- dplyr::case_when(
    d$alpha_group == "耗竭型"  &  d$self_ok ~ "潜在抑郁",
    d$alpha_group == "耗竭型"  & !d$self_ok ~ "典型抑郁",
    d$alpha_group == "高意识型" &  d$self_ok ~ "真正健康",
    TRUE ~ NA_character_
  )
  d$potential_depression <- d$self_ok & d$alpha_group == "耗竭型"   # 论文「潜在痛苦」
  d$hidden_distress      <- d$potential_depression                  # 旧代码同义名
  d$healthy_control      <- d$self_ok & d$alpha_group == "高意识型"
  d
}

# 书籍已公开的唯一下游 HCF 变换：H_trend = α1²。
# 不是额外后台读数或实测纵向变化；不派生其他指标，不重构后台 H。
add_h_terms <- function(d) {
  d$H_trend <- d$alpha1^2
  d
}

# ---------- BMI：单位修正 + WHO 阈值 ----------
BMI_OBESE_CUT <- 30   # 实证：阈值 30 与 Prolific obese 标签 471/471 一致
bmi_from_hw <- function(height_raw, weight_raw) {
  h <- suppressWarnings(as.numeric(height_raw))
  w <- suppressWarnings(as.numeric(weight_raw))
  h <- ifelse(!is.na(h) & h < 3, h * 100, h)          # 有人填 1.64（米）
  h <- ifelse(!is.na(h) & (h < 100 | h > 250), NA_real_, h)
  w <- ifelse(!is.na(w) & (w < 30 | w > 250), NA_real_, w)
  b <- ifelse(!is.na(h) & !is.na(w), w / (h / 100)^2, NA_real_)
  ifelse(!is.na(b) & (b < 10 | b > 80), NA_real_, b)
}
