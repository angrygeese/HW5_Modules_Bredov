import pickle
import random


import multiprocessing
import numpy as np
import os
import pandas as pd
import threading


from concurrent.futures import (ProcessPoolExecutor, ThreadPoolExecutor)
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier


SEED = 111
random.seed(SEED)
np.random.seed(SEED)


class RandomForestClassifierCustomMT(BaseEstimator):
    def __init__(
        self, n_estimators=10,
        max_depth=None,
        max_features=None,
        random_state=SEED,
        tmp_dir='example_data/processes_dumps'
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.tmp_dir = tmp_dir  # https://stackoverflow.com/a/45024475

        self.trees = []
        self.feat_ids_by_tree = []
        self.probas_y = []

    @property
    def tmp_dir(self):
        return self._tmp_dir

    @tmp_dir.setter
    def tmp_dir(self, tmp_dir):
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        self._tmp_dir = tmp_dir

    def fit_dts(self, X, y, estimators):
        n_samples, n_features = X.shape

        for tree_num in estimators:
            np.random.seed(seed=self.random_state + tree_num)
            feature_idx = np.random.choice(
                n_features,
                self.max_features,
                replace=False
                )

            btstrp_sample = np.random.choice(
                n_samples,
                size=n_samples,
                replace=True
                )
            btstrp_sample_X = X[btstrp_sample][:,feature_idx]
            btstrp_sample_y = y[btstrp_sample]

            model = DecisionTreeClassifier(
                max_depth=self.max_depth,
                max_features=self.max_features,
                random_state=self.random_state

                )
            model.fit(btstrp_sample_X, btstrp_sample_y)

            with open(f'{self.tmp_dir}/model_{tree_num}_dump.pkl', 'wb') as md:
                pickle.dump((model, feature_idx), md)  # https://stackoverflow.com/a/19201448

        return self

    def split_estimators(self, n_jobs):
        floordiv, mod = self.n_estimators // n_jobs, self.n_estimators % n_jobs
        idx = np.arange(floordiv, floordiv * n_jobs, floordiv)
        if mod:
            idx += np.append(
                    np.zeros(n_jobs - mod, dtype='int8'),
                    np.arange(1, mod))
        return np.split(np.arange(self.n_estimators), idx)

    def assemble_data_from_process(self, estimators, name='model'):
        for tree_num in estimators:
            with open(f'{self.tmp_dir}/{name}_{tree_num}_dump.pkl', 'rb') as md:
                match name:
                    case 'model':
                        model, model_feature_idx = pickle.load(md)
                        self.trees.append(model)
                        self.feat_ids_by_tree.append(model_feature_idx)
                    case 'proba':
                        self.probas_y.append(pickle.load(md))

    def fit(self, X, y, n_jobs):
        self.classes_ = sorted(np.unique(y))

        processes = [multiprocessing.Process(
            target=self.fit_dts, args=(X, y, estimators)
            ) for estimators in self.split_estimators(n_jobs)]
        for proc in processes:
            proc.start()
        for proc in processes:
            proc.join()

        threads = [threading.Thread(
            target=self.assemble_data_from_process, args=(estimators, 'model')
            ) for estimators in self.split_estimators(n_jobs)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        return self

    def predict_proba_dts(self, X, estimators):

        for tree_num in estimators:
            tree = self.trees[tree_num]
            tree_proba = tree.predict_proba(X[:, self.feat_ids_by_tree[tree_num]])

            with open(f'{self.tmp_dir}/proba_{tree_num}_dump.pkl', 'wb') as md:
                pickle.dump(tree_proba, md)

    def predict_proba(self, X, n_jobs):
        self.probas_y.clear()
        
        threads = [threading.Thread(
            target=self.predict_proba_dts, args=(X, estimators)
            ) for estimators in self.split_estimators(n_jobs)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        threads = [threading.Thread(
            target=self.assemble_data_from_process, args=(estimators, 'proba')
            ) for estimators in self.split_estimators(n_jobs)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        return np.array(self.probas_y).mean(axis=0)

    def predict(self, X, n_jobs):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)
        self.probas_y.clear()

        return predictions


class RandomForestClassifierCustomST(BaseEstimator):
    def __init__(
        self, n_estimators: int = 10, max_depth=None, max_features=None, random_state=SEED
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit_tree(self, X, y, estimator):
        n_samples, n_features = X.shape

        np.random.seed(seed=self.random_state + estimator)
        feature_idx = np.random.choice(
            n_features,
            self.max_features,
            replace=False
            )

        btstrp_sample = np.random.choice(
            n_samples,
            size=n_samples,
            replace=True
            )
        btstrp_sample_X = X[btstrp_sample][:,feature_idx]
        btstrp_sample_y = y[btstrp_sample]

        model = DecisionTreeClassifier(
            max_depth=self.max_depth,
            max_features=self.max_features,
            random_state=self.random_state
            )
        model.fit(btstrp_sample_X, btstrp_sample_y)

        return model, feature_idx

    def fit(self, X, y, n_jobs):
        self.classes_ = sorted(np.unique(y))

        args = [(X, y, estimator) for estimator in range(self.n_estimators)]
        with ProcessPoolExecutor(n_jobs) as pool:
            self.trees, self.feat_ids_by_tree = zip(*list(pool.map(self.fit_tree, *zip(*args))))  # https://stackoverflow.com/a/6976772
                                                                                                  # https://stackoverflow.com/a/62440230
        return self

    def predict_proba_dt(self, X, estimator):
        tree = self.trees[estimator]
        tree_proba = tree.predict_proba(X[:, self.feat_ids_by_tree[estimator]])

        return tree_proba

    def predict_proba(self, X, n_jobs):
        args = [(X, estimator) for estimator in range(self.n_estimators)]
        with ThreadPoolExecutor(n_jobs) as pool:
            probas_y = list(pool.map(self.predict_proba_dt, *zip(*args)))

        return np.array(probas_y).mean(axis=0)

    def predict(self, X, n_jobs):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
