# LDPC-BIKE

Replace BIKE's decoder with sum-product and min-sum.

---

在 LEVEL-1，sum-product 表现较好，dfr 接近 0，但是效率比 BGF 低很多。

当 error 的 weight 小于 100 时，min-sum decoder 表现较好，但是在 LEVEL-1 (error 的 weight 为 134) 表现很差。