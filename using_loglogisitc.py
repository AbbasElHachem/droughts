
#     def fisk_of_pet_data(self, pet_data):
#         c, loc, scale = stats.fisk.fit(pet_data)
#         return stats.fisk.cdf(pet_data, c=c, scale=scale,
#                               loc=loc)
#
# #     def likelihood_log_logistic(self, params, _data):
# #         return -np.sum(np.vectorize(self.fisk_of_pet_data)
# #                        (_data, params[0], params[1], params[2]))
#
# #     def optimize_likelihood_log_logistic(self, _data):
# #
# #         _params = optimize.minimize(self.likelihood_log_logistic,
# #                                     args=(_data),
# #                                     x0=[1, 0.1, 2], method='SLSQP',
# #                                     bounds=[(0., None), (-1., 1.), (None, None)])
# #         return _params
#
#     def log_logistic_prob(self, vals, *params):
#         return(1 + (params[0] / (vals - params[2]))**params[1])**-1
#
#     def spei_calculation(self, cdf_vals):
#         c0, c1, c2 = 2.515517, 0.802853, 0.010328
#         d1, d2, d3 = 1.432788, 0.189269, 0.001308
#         SPEI_total = np.empty(shape=(len(cdf_vals)))
#
#         for i, val in enumerate(cdf_vals):
#             try:
#                 if 0 < val <= 0.5:
#                     w = -2 * np.log(val)
#                     SPEI_total[i] = -(w - (c0 + c1 * w + c2 * w**2) /
#                                       (1 + d1 * w + d2 * w**2 + d3 * w**3))
#                 elif 0.5 <= val < 1:
#                     w = -2 * np.log(1 - val)
#
#                     SPEI_total[i] = (w - (c0 + c1 * w + c2 * w**2) /
#                                      (1 + d1 * w + d2 * w**2 + d3 * w**3))
#
#             except Exception:
#                 SPEI_total[i] = np.nan

#         return SPEI_total
