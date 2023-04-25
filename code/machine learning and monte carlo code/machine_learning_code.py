import pandas as pd
df = pd.read_csv("MgFe2O4-dataset.csv")
X = df.iloc[:, 1:24]
y = df.iloc[:, 24]/8

from sklearn.pipeline import make_pipeline
#import train_test_split library
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.compose import make_column_transformer
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
import numpy as np
import matplotlib.pyplot as plt 
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_validate
import pickle
from joblib import dump, load


"""## Linear Regression"""

Model = linear_model.LinearRegression(fit_intercept=True)
cv = ShuffleSplit(n_splits=4, test_size=0.25)
scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=4, return_train_score=True)
test_RMSE = -np.average(scores['test_score'])
train_RMSE = -np.average(scores['train_score'])
dump(Model, 'lr-MgFe.joblib') 
print(test_RMSE)
print(train_RMSE)



"""## Ridge Regression"""

from sklearn.linear_model import Ridge

Model = Ridge(alpha=21.5)
scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=4, return_train_score=True)
test_RMSE = -np.average(scores['test_score'])
train_RMSE = -np.average(scores['train_score'])
print(test_RMSE)
print(train_RMSE)


"""# Support Vector Machine"""

from sklearn.svm import SVR
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

Model = make_pipeline(SVR(C=2.3, epsilon=0, kernel = 'rbf'))
scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=cv, return_train_score=True)
test_RMSE = -np.average(scores['test_score'])
train_RMSE = -np.average(scores['train_score'])
dump(Model, 'svm-MgFe.joblib') 
print(test_RMSE, train_RMSE)


for c in np.arange(1, 3, 0.1):
  Model = make_pipeline(SVR(C=c, epsilon=0, kernel = 'sigmoid'))
  cv = ShuffleSplit(n_splits=4, test_size=0.25)
  scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=cv, return_train_score=True)
  test_RMSE = -np.average(scores['test_score'])
  train_RMSE = -np.average(scores['train_score'])
  # print(test_RMSE, train_RMSE)
  print("C = {}, epsilon = {}, test error = {}, train error = {}".format(c, 0, test_RMSE, train_RMSE))


"""# Neural Network"""

from sklearn.neural_network import MLPRegressor

Model = make_pipeline(MLPRegressor(learning_rate = 'adaptive', learning_rate_init = 0.01, max_iter=10000, activation = "logistic", solver = "lbfgs", alpha = 0.01, hidden_layer_sizes = (120, 160)))
cv = ShuffleSplit(n_splits=4, test_size=0.25)
scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=cv, return_train_score=True)
test_RMSE = -np.average(scores['test_score'])
train_RMSE = -np.average(scores['train_score'])
dump(Model, 'nn-MgFe.joblib')

"""# Decision tree"""

from sklearn.tree import DecisionTreeRegressor

Model = make_pipeline(DecisionTreeRegressor(criterion='absolute_error', splitter = "best", max_features = 10, max_depth=5))
cv = ShuffleSplit(n_splits=4, test_size=0.25)
scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=4, return_train_score=True)
test_RMSE = -np.average(scores['test_score'])
train_RMSE = -np.average(scores['train_score'])
print(test_RMSE, train_RMSE)

print(scores, np.average(scores))

"""## Random Forest"""

from sklearn.ensemble import RandomForestRegressor

Model = make_pipeline(RandomForestRegressor(n_estimators=2500,criterion='absolute_error',max_depth = 8 ,max_features = 4))
cv = ShuffleSplit(n_splits=4, test_size=0.25)
scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=cv, return_train_score=True)
test_RMSE = -np.average(scores['test_score'])
train_RMSE = -np.average(scores['train_score'])
print(test_RMSE, train_RMSE)

"""# Gradient Boosting Machine"""

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.inspection import permutation_importance

cv = ShuffleSplit(n_splits=4, test_size=0.25)
scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=cv, return_train_score=True)
test_RMSE = -np.average(scores['test_score'])
train_RMSE = -np.averagModel = make_pipeline(GradientBoostingRegressor(subsample = 0.5, max_depth = 7, min_samples_leaf =5, min_samples_split = 4, loss = 'squared_error', validation_fraction = 0.9, learning_rate = 0.01, alpha = 0.01, n_estimators = 7000))
e(scores['train_score'])
dump(Model, 'gbm-MgFe.joblib') 
print(test_RMSE, train_RMSE)


Model = GradientBoostingRegressor(subsample = 0.5, max_depth = 7, min_samples_leaf =5, min_samples_split = 4, loss = 'squared_error', validation_fraction = 0.9, learning_rate = 0.01, alpha = 0.01, n_estimators = 7000)
feature_importance = Model.fit(X, y).feature_importances_
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + 0.5
fig = plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.barh(pos, feature_importance[sorted_idx], align="center")
plt.yticks(pos, np.array(X.columns)[sorted_idx])
plt.title("Feature Importance (MDI)")



"""## K-Nearest Neighbor"""

from sklearn import neighbors

for k in range(1, 20):
  Model = neighbors.KNeighborsRegressor(n_neighbors = k)
  cv = ShuffleSplit(n_splits=4, test_size=0.25)
  scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=4, return_train_score=True)
  test_RMSE = -np.average(scores['test_score'])
  train_RMSE = -np.average(scores['train_score'])
  print(test_RMSE, train_RMSE)

"""# Error with size"""
train_error = []
test_error = []

for i in range(100):
  df_new = df.sample(n = 10)
  X = df_new.iloc[:,1:24]
  y = df_new.iloc[:,24]/8
  Model = make_pipeline(SVR(C=2.3, epsilon=0, kernel='rbf'))
  cv = ShuffleSplit(n_splits=4, test_size=0.25)
  scores = cross_validate(Model, X, y, scoring='neg_root_mean_squared_error', cv=cv, return_train_score=True)
  test_RMSE = -np.average(scores['test_score'])
  train_RMSE = -np.average(scores['train_score'])
  train_error.append(train_RMSE)
  test_error.append(test_RMSE)
print(np.average(train_error), np.average(test_error))

