import cmocean
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage.filters import gaussian_filter
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from read_scenario import *
import tikzplotlib
import tomli
import scipy

matplotlib.use('Qt5Agg')
mumblue = np.array([0, 100, 222]) / 255
mumred = np.array([220, 33, 77]) / 255
mumred_alt = np.array([220, 33, 77, 155]) / 255
mumgreen = np.array([0, 140, 0]) / 255


@dataclass
class EstimatePlotConfig():
    plot_trajectory: bool = True
    plot_agent: bool = True
    plot_points: bool = True


class GTInterface(object):
    def __init__(self, name):
        self.path = os.path.join("../../resources_build/configs/maps/", name)
        file = open(self.path, "rb")
        toml_dict = tomli.load(file)
        self.sx = toml_dict["general"]["size_x"]
        self.sy = toml_dict["general"]["size_y"]
        self.values = np.array(toml_dict["value_field"]["values"])
        self.dx = self.sx / (self.values.shape[0] - 1)
        self.dy = self.sy / (self.values.shape[1] - 1)

    def get_value(self, x, y):
        # interpolate value based on size etc., first find closest grid point
        id_x = x / self.dx
        id_y = y / self.dy
        # linear interpolation
        x1 = int(np.floor(id_x))
        x2 = x1 + 1
        y1 = int(id_y)
        y2 = y1 + 1
        if x2 >= self.values.shape[0]:
            x2 = x1
        if y2 >= self.values.shape[1]:
            y2 = y1
        v11 = self.values[x1, y1]
        v12 = self.values[x1, y2]
        v21 = self.values[x2, y1]
        v22 = self.values[x2, y2]
        w11 = (x2 * self.dx - x) * (y2 * self.dy - y)
        w12 = (x2 * self.dx - x) * (y - y1 * self.dy)
        w21 = (x - x1 * self.dx) * (y2 * self.dy - y)
        w22 = (x - x1 * self.dx) * (y - y1 * self.dy)
        v = (v11 * w11 + v12 * w12 + v21 * w21 + v22 * w22) / (
                w11 + w12 + w21 + w22)
        return v


class SingleExperiment(object):
    def __init__(self, path):
        self.metrics = None
        self.metrics_keys = None
        self.path = path
        self.num_steps = 0

        config_file = open(os.path.join(self.path, "config.txt"), "r")
        self.field_name = config_file.readline()
        print("Reading " + self.path + " with Scenario " + self.field_name)
        config_file.close()
        self.getTimeSteps()
        self.read_metrics()
        gmrf_conf = open(os.path.join(self.path, "configs", "gmrf_conf.toml"),
                         "rb")
        toml_dict = tomli.load(gmrf_conf)
        self.nx = toml_dict["general"]["N_X"]
        self.ny = toml_dict["general"]["N_Y"]
        self.estimates = {}
        self.covariances = {}
        self.measurement_locs = None
        # self.read_data()

    def get_source_map(self, t_id):
        source_file = open(os.path.join(self.path, "source_estimates.csv"), "r")
        source_csv = csv.reader(source_file, delimiter=",")
        source = np.zeros((self.nx, self.ny))
        c = -3
        if t_id == -1:
            t_id = self.num_steps
        for source_read in source_csv:
            c += 1
            if c == t_id - 1:
                all_sources = np.array(source_read[1:-1], dtype=np.float64)
                l = len(all_sources)
                return np.reshape(all_sources, (self.ny, self.nx))

    def getErrorMap(self, t):
        gt = ScenarioLoader(
            "../../resources_build/configs/maps/" + self.field_name + ".toml")
        gt.load_scenario()
        est_file = open(os.path.join(self.path, "gmrf_estimates.csv"), "r")
        est_csv = csv.reader(est_file, delimiter=",")
        est = np.zeros((self.nx, self.ny))
        c = -3
        if t == -1:
            t = self.num_steps
        for est_read in est_csv:
            c += 1
            if c == t - 1:
                all_ests = np.array(est_read[1:-1], dtype=np.float64)
                l = len(all_ests)
                k = round(l / 3)
                est = np.clip(
                    np.reshape(all_ests[:k], (self.ny, self.nx)).astype(
                        np.float64).T, 0, 1)
                break
        est_file.close()
        vals = np.zeros((self.nx, self.ny))
        for i in range(self.nx):
            for j in range(self.ny):
                vals[i, j] = gt.get_value(i / self.nx * gt.sx,
                                          j / self.ny * gt.sy)
        return np.abs(est - vals)

    def getTimeSteps(self):
        est_file = open(os.path.join(self.path, "gmrf_estimates.csv"), "r")
        num_lines = sum(1 for _ in est_file)
        self.num_steps = num_lines - 2
        est_file.close()

    def read_data(self):
        self.estimates["field"] = np.zeros((self.nx, self.ny, self.num_steps))
        self.estimates["flow_x"] = np.zeros((self.nx, self.ny, self.num_steps))
        self.estimates["flow_y"] = np.zeros((self.nx, self.ny, self.num_steps))
        self.covariances["field"] = np.zeros((self.nx, self.ny, self.num_steps))
        self.covariances["flow_x"] = np.zeros(
            (self.nx, self.ny, self.num_steps))
        self.covariances["flow_y"] = np.zeros(
            (self.nx, self.ny, self.num_steps))

        est_file = open(os.path.join(self.path, "gmrf_estimates.csv"), "r")
        est_csv = csv.reader(est_file, delimiter=",")

        cov_file = open(os.path.join(self.path, "field_covariances.csv"), "r")
        cov_csv = csv.reader(cov_file, delimiter=",")

        flow_x_cov_file = open(os.path.join(self.path, "flowx_covariances.csv"),
                               "r")
        flow_x_cov_csv = csv.reader(flow_x_cov_file, delimiter=",")

        flow_y_cov_file = open(os.path.join(self.path, "flowy_covariances.csv"),
                               "r")
        flow_y_cov_csv = csv.reader(flow_y_cov_file, delimiter=",")
        c = -3
        for est, cov, flow_x_cov, flow_y_cov in zip(est_csv, cov_csv,
                                                    flow_x_cov_csv,
                                                    flow_y_cov_csv):
            c += 1
            if c >= 0:
                all_ests = est[1: -1]
                l = len(all_ests)
                k = round(l / 3)

                covs = cov[1: -1]
                flow_x_covs = flow_x_cov[1: -1]
                flow_y_covs = flow_y_cov[1: -1]

                mean_ests = np.clip(
                    np.reshape(all_ests[:k], (self.ny, self.nx)).astype(
                        np.float64).T, 0, 1)
                flow_x_ests = np.reshape(all_ests[k: 2 * k],
                                         (self.ny, self.nx)).astype(
                    np.float64).T
                flow_y_ests = np.reshape(all_ests[2 * k: 3 * k],
                                         (self.ny, self.nx)).astype(
                    np.float64).T
                self.estimates["field"][:, :, c] = mean_ests
                self.estimates["flow_x"][:, :, c] = flow_x_ests
                self.estimates["flow_y"][:, :, c] = flow_y_ests

                covs = np.reshape(covs, (self.ny, self.nx)).astype(np.float64).T
                flow_x_covs = np.reshape(flow_x_covs,
                                         (self.ny, self.nx)).astype(
                    np.float64).T
                flow_y_covs = np.reshape(flow_y_covs,
                                         (self.ny, self.nx)).astype(
                    np.float64).T
                self.covariances["field"][:, :, c] = covs
                self.covariances["flow_x"][:, :, c] = flow_x_covs
                self.covariances["flow_y"][:, :, c] = flow_y_covs

        print("Done reading!")

    def read_metrics(self):
        self.metrics = {}
        df = pd.read_csv(os.path.join(self.path, "metrics.csv"),
                         engine="pyarrow")
        for a in df.columns:
            self.metrics[a] = df[a].values
        self.metrics["FieldLogLikelihoodMetric"] = -self.metrics[
            "FieldLogLikelihoodMetric"]
        # self.metrics["source_threshold"] = [np.argmax(
        #     self.metrics["SourceDistMetric"] < 1.0)]
        self.metrics_keys = list(self.metrics.keys())

    def plot_source(self, axs):
        source_file = open(os.path.join(self.path, "Source.csv"), "r")
        csv_reader = csv.reader(source_file, delimiter=",")
        c = 0
        x_vals = []
        y_vals = []
        for row in csv_reader:
            if c == self.num_steps:
                x_vals.append(float(row[1]))
                y_vals.append(float(row[2]))
            c += 1
        axs.scatter(x_vals, y_vals, c="white", s=100, marker="+")

    def get_sources(self, t=-1):
        source_file = open(os.path.join(self.path, "Source.csv"), "r")
        csv_reader = csv.reader(source_file, delimiter=",")
        c = 0
        x_vals = []
        y_vals = []
        t_range = []
        if t == "All":
            t_range = [0, self.metrics["t"][-1]]
        elif isinstance(t, list):
            t_range = [t[0] + 1, t[1] + 1]
        elif t == -1:
            t_range = [self.metrics["t"][-1], self.metrics["t"][-1]]
        else:
            t_range = [t, t]
        for row in csv_reader:
            if c > 1:
                if t_range[0] <= self.metrics["t"][c - 1] <= t_range[1]:
                    x_vals.append(float(row[1]))
                    y_vals.append(float(row[2]))
            c += 1
        return x_vals, y_vals

    def get_velocities_with_locations(self, ):
        states_file = open(os.path.join(self.path, "state.csv"), "r")
        csv_reader = csv.reader(states_file, delimiter=",")
        c = 0
        x_vals = []
        y_vals = []
        v_norm_vals = []
        for row in csv_reader:
            if c > 1:
                num_a = int((len(row) - 1) / 4)
                for i in range(num_a):
                    x_vals.append(float(row[1 + 4 * i]))
                    y_vals.append(float(row[2 + 4 * i]))
                    v_norm_vals.append((float(row[4 + 4 * i]))**2)
            c += 1
        x_vals = np.array(x_vals)
        # reshape to (num_agents, num_steps)
        x_vals = np.reshape(x_vals, (self.num_steps,
                                     int(len(x_vals) / self.num_steps))).T
        y_vals = np.array(y_vals)
        y_vals = np.reshape(y_vals, (self.num_steps,
                                     int(len(y_vals) / self.num_steps))).T
        v_norm_vals = np.array(v_norm_vals)
        v_norm_vals = np.reshape(v_norm_vals, (self.num_steps,
                                                  int(len(v_norm_vals) / self.num_steps))).T
        return x_vals, y_vals, v_norm_vals

    def get_states(self, rw=False):
        states_file = open(os.path.join(self.path, "state.csv"), "r")
        csv_reader = csv.reader(states_file, delimiter=",")
        c = 0
        x_vals = []
        y_vals = []
        for row in csv_reader:
            if c >= 2:
                if rw:
                    num_a = int((len(row) - 1) / 3)
                    for i in range(num_a):
                        x_vals.append(float(row[1 + 3 * i]))
                        y_vals.append(float(row[2 + 3 * i]))
                else:
                    if len(row) == 3:
                        x_vals.append(float(row[1]))
                        y_vals.append(float(row[2]))
                    else:
                        num_a = int((len(row) - 1) / 4)
                        for i in range(num_a):
                            x_vals.append(float(row[1 + 4 * i]))
                            y_vals.append(float(row[2 + 4 * i]))
            c += 1
        x_vals = np.array(x_vals)
        # reshape to (num_agents, num_steps)
        x_vals = np.reshape(x_vals, (self.num_steps,
                                     int(len(x_vals) / self.num_steps))).T
        y_vals = np.array(y_vals)
        y_vals = np.reshape(y_vals, (self.num_steps,
                                     int(len(y_vals) / self.num_steps))).T

        return x_vals, y_vals

    def plot_timestep(self, step_id, export=None, axs=None, config=None):
        plot_agent = False
        plot_trajectory = True
        plot_points = True
        if config is not None:
            plot_agent = config.plot_agent
            plot_trajectory = config.plot_trajectory
            plot_points = config.plot_points
        x, y = self.get_states()
        scen = ScenarioLoader(
            "../../resources_build/configs/maps/" + self.field_name + ".toml")
        scen.load_scenario()
        selfplot = False
        if axs is None:
            fig, axs = plt.subplots(1, 1, figsize=(10, 10))
            selfplot = True
        scen.plot_value_field(axs, False, False)
        sources = self.get_sources(step_id)
        axs.set_title("t = " + str(self.metrics["t"][step_id]))

        if plot_agent:
            cirle = plt.Circle((x[step_id], y[step_id]), 0.5, color="black")
            axs.add_artist(cirle)
        axs.scatter(sources[0], sources[1], c="black", s=1000, marker="*")
        axs.scatter(sources[0], sources[1], c="white", s=400, marker="*")
        axs.set_aspect("equal")
        if plot_trajectory:
            axs.plot(x[:step_id], y[:step_id], c="black")
        if plot_points:
            colors = cm.plasma(np.linspace(0, 1, len(x[:step_id])))
            axs.scatter(x[:step_id], y[:step_id], c=colors)
        if export is not None:
            plt.savefig(str(export) + ".svg")
        if selfplot:
            plt.show()

    def render_video(self):
        for i in range(self.num_steps):
            plot_conf = EstimatePlotConfig()
            plot_conf.plot_agent = True
            plot_conf.plot_trajectory = True
            plot_conf.plot_points = False
            fig, axs = plt.subplots(1, 1, figsize=(20, 10))
            # self.plot_timestep(i, axs=axs[0], config=plot_conf)
            self.plot_estimates(axs, i)
            self.plot_timestep(i, axs=axs, config=plot_conf)
            # add the leading zeros to the name, could be thousands of pictures
            name = str(i)
            while len(name) < len(str(self.num_steps)):
                name = "0" + name
            plt.savefig(os.path.join("video", name + ".png"))
            plt.close(fig)

    def plot_mean_error(self, axs, t=-1):
        gt = ScenarioLoader(
            "../../resources_build/configs/maps/" + self.field_name + ".toml")
        gt.load_scenario()
        vals = np.zeros((self.nx, self.ny))
        for i in range(self.nx):
            for j in range(self.ny):
                vals[i, j] = gt.get_value(i / self.nx * gt.sx,
                                          j / self.ny * gt.sy)
        axs.imshow(np.abs(vals - self.estimates["field"][:, :, t]).T,
                   cmap="plasma", origin="lower",
                   extent=[0, gt.sx, 0, gt.sy], vmin=0, vmax=1,
                   interpolation="bicubic")
        print(np.max(np.abs(vals - self.estimates["field"][:, :, t])))
        gt.plot_walls(axs)

    def plot_evaluation(self, axs, t=-1):
        gt = ScenarioLoader(
            "../../resources_build/configs/maps/" + self.field_name + ".toml")
        gt.load_scenario()
        self.plot_estimates(axs[0, 0], t)
        axs[0, 0].set_title("t = " + str(self.metrics["t"][t]))
        axs[0, 0].set_aspect("equal")
        # quiver with flow
        # create meshgrid
        x = np.linspace(0, gt.scenario.general.size_x, self.nx)
        y = np.linspace(0, gt.scenario.general.size_y, self.ny)
        X, Y = np.meshgrid(x, y)
        x_ = np.linspace(0, gt.scenario.general.size_x, 50)
        y_ = np.linspace(0, gt.scenario.general.size_y, 25)
        X_, Y_ = np.meshgrid(x_, y_)
        u_interp = interpolate.griddata((X.flatten(), Y.flatten()),
                                        self.estimates["flow_x"][:, :,
                                        t].T.flatten(), (X_, Y_))
        v_interp = interpolate.griddata((X.flatten(), Y.flatten()),
                                        self.estimates["flow_y"][:, :,
                                        t].T.flatten(), (X_, Y_))
        norm = np.sqrt(u_interp ** 2 + v_interp ** 2)
        u_interp = u_interp / np.where(norm == 0, 1, norm)
        v_interp = v_interp / np.where(norm == 0, 1, norm)
        # plot quiver
        axs[1, 0].quiver(X_, Y_, u_interp, v_interp, color="black", scale=30)
        # axs[1].scatter(measurement_locs[:, 0], measurement_locs[:, 1], c="black", s=100, marker="+")
        axs[1, 0].set_aspect("equal")
        axs[1, 0].set_title("flow")
        axs[0, 1].imshow(self.covariances["field"][:, :, t].T, cmap="plasma",
                         origin="lower", extent=[0, gt.sx, 0, gt.sy])
        axs[0, 1].set_aspect("equal")
        axs[0, 1].set_title("covariances")
        axs[1, 1].imshow((self.covariances["flow_x"][:, :, t] +
                          self.covariances["flow_y"][:, :, t]).T / 2,
                         cmap="plasma", origin="lower",
                         extent=[0, gt.sx, 0, gt.sy])
        axs[1, 1].set_aspect("equal")
        axs[1, 1].set_title("covariances")
        likelihoods = np.zeros((self.nx, self.ny))
        vals = np.zeros((self.nx, self.ny))
        log_likelihoods = np.zeros((self.nx, self.ny))
        for i in range(self.nx):
            for j in range(self.ny):
                likelihoods[i, j] = 1 / (np.sqrt(
                    self.covariances["field"][i, j, t] * 2 * np.pi)) * np.exp(
                    -0.5 * (gt.get_value(i / self.nx * gt.sx,
                                         j / self.ny * gt.sy) -
                            self.estimates["field"][i, j, t]) ** 2 /
                    self.covariances["field"][i, j, t])
                log_likelihoods[i, j] = - (gt.get_value(i / self.nx * gt.sx,
                                                        j / self.ny * gt.sy) -
                                           self.estimates["field"][
                                               i, j, t]) ** 2 / \
                                        self.covariances["field"][
                                            i, j, t] - np.log(
                    np.sqrt(self.covariances["field"][i, j, t] * 2 * np.pi))
                vals[i, j] = gt.get_value(i / self.nx * gt.sx,
                                          j / self.ny * gt.sy)
        axs[2, 0].imshow(np.abs(vals - self.estimates["field"][:, :, t]).T,
                         cmap="viridis", origin="lower",
                         extent=[0, gt.sx, 0, gt.sy], vmin=0, vmax=1)
        axs[2, 0].set_title("Mean Error")
        print("Mean Error :",
              np.sqrt(np.sum((vals - self.estimates["field"][:, :,
                                     t]).T ** 2) / self.nx / self.ny))
        print("Mean Likelihood :", np.mean(likelihoods))
        print("Joint Likelihood :", np.sum(log_likelihoods))
        print("Min Likelihood :", np.min(log_likelihoods))
        print("Max Likelihood :", np.max(log_likelihoods))
        axs[2, 1].imshow(log_likelihoods.T, cmap="plasma", origin="lower",
                         extent=[0, gt.sx, 0, gt.sy])
        axs[2, 1].set_aspect("equal")
        axs[2, 1].set_title("Log Likelihood")

        gt.plot_walls(axs[0, 0])
        gt.plot_walls(axs[0, 1])
        gt.plot_walls(axs[1, 0])
        gt.plot_walls(axs[1, 1])
        gt.plot_walls(axs[2, 0])
        gt.plot_walls(axs[2, 1])

    def plot_estimates(self, axs, t=-1):
        gt = ScenarioLoader(
            "../../resources_build/configs/maps/" + self.field_name + ".toml")
        gt.load_scenario()

        axs.imshow(self.estimates["field"][:, :, t].T, cmap="plasma",
                   origin="lower", extent=[0, gt.scenario.general.size_x, 0,
                                           gt.scenario.general.size_y],
                   interpolation="bicubic")


class RepeatedExperiment(object):
    def __init__(self, path, name):
        self.name = name
        self.path = path
        self.experiments = []
        self.field_name = None
        self.mean_metrics = None
        self.metrics_keys = None
        self.varying_lengths = False
        for directory in os.listdir(self.path):
            self.experiments.append(
                SingleExperiment(os.path.join(self.path, directory)))
            if self.field_name is None:
                self.field_name = self.experiments[-1].field_name
            elif not self.field_name == self.experiments[-1].field_name:
                raise RuntimeError("Eval Fields do not match up.")
        self.mean_metrics = {}
        self.median_metrics = {}
        self.std_metrics = {}
        self.q1_metrics = {}
        self.q2_metrics = {}
        self.metrics = {}
        self.t = None
        self.t_ = []
        self.get_metrics()

    def get_final_error(self):
        source_dists = self.metrics["SourceDistMetric"]
        times = []
        errors = []
        for i, v in enumerate(source_dists):
            times.append(len(v))
            errors.append(v[-1])
        times = np.array(times)
        errors = np.array(errors)
        return errors, times

    def plot_time_visited_heatmap(self, axs, export=None, t_max=-1):
        x = []
        y = []
        t = []
        map_name = self.experiments[0].field_name
        example_x = None
        example_y = None
        example_sources_x = None
        example_sources_y = None
        binsx = self.experiments[0].nx * 2
        binsy = self.experiments[0].ny * 2
        scen = ScenarioLoader(
            "../../resources_build/configs/maps/" + map_name + ".toml")
        scen.load_scenario()
        for i, exp in enumerate(self.experiments):
            t_vals = exp.metrics["t"]
            if t_max == -1:
                t_max = t_vals[-1]
            t_id = np.argmin(np.abs(t_vals - t_max))
            x_v, y_v, v_v = exp.get_velocities_with_locations()
            x_v = x_v[:, :t_id]
            y_v = y_v[:, :t_id]
            t_ids = np.arange(0, t_id) # times number of agents into array like above
            # reverse t_ids
            t_ids = t_ids[::-1]
            t_ids = np.repeat(t_ids[np.newaxis, :], x_v.shape[0], axis=0)
            x = x + x_v.flatten().tolist()
            y = y + y_v.flatten().tolist()
            t = t + t_ids.flatten().tolist()

        x = np.array(x)
        y = np.array(y)
        t = np.array(t)
        # else:
        hist, xedges, yedges = np.histogram2d(x, y,
                                              bins=[binsx,
                                                    binsy],
                                              range=[[0, scen.sx],
                                                     [0, scen.sy]],
                                              weights=t)
        count, _, _ = np.histogram2d(x, y, bins=[binsx,
                                                 binsy],
                                     range=[[0, scen.sx],
                                            [0, scen.sy]])
        with np.errstate(divide='ignore', invalid='ignore'):
            average_speeds = np.true_divide(hist, count)
            average_speeds[~np.isfinite(average_speeds)] = 0
        # heatmap = gaussian_filter(heatmap, sigma=1.2)

        extent = [0, scen.sx, 0, scen.sy]
        # plt gt
        print(average_speeds.sum())
        scen.plot_value_field(axs, False, False)
        axs.imshow(average_speeds.T, extent=extent, origin="lower",
                   cmap="viridis", vmin=0, vmax=t_id)
        if export is not None:
            plt.savefig(export + ".svg", dpi=400)
        # plt.show()

    def plot_velocity_heatmap(self, axs, export=None, t_max=-1):
        x = []
        y = []
        v = []
        map_name = self.experiments[0].field_name
        example_x = None
        example_y = None
        example_sources_x = None
        example_sources_y = None
        binsx = self.experiments[0].nx * 2
        binsy = self.experiments[0].ny * 2
        scen = ScenarioLoader(
            "../../resources_build/configs/maps/" + map_name + ".toml")
        scen.load_scenario()
        for i, exp in enumerate(self.experiments):
            t_vals = exp.metrics["t"]
            if t_max == -1:
                t_max = t_vals[-1]
            t_id = np.argmin(np.abs(t_vals - t_max))
            x_v, y_v, v_v = exp.get_velocities_with_locations()
            x_v = x_v[:, :t_id]
            y_v = y_v[:, :t_id]
            v_v = v_v[:, :t_id]

            x = x + x_v.flatten().tolist()
            y = y + y_v.flatten().tolist()
            v = v + v_v.flatten().tolist()

        x = np.array(x)
        y = np.array(y)
        # else:
        hist, xedges, yedges = np.histogram2d(x, y,
                                                 bins=[binsx,
                                                       binsy],
                                                 range=[[0, scen.sx],
                                                        [0, scen.sy]],
                                                 weights=v)
        count, _, _ = np.histogram2d(x, y, bins=[binsx,
                                                    binsy],
                                        range=[[0, scen.sx],
                                               [0, scen.sy]])
        with np.errstate(divide='ignore', invalid='ignore'):
            average_speeds = np.true_divide(hist, count)
            average_speeds[~np.isfinite(average_speeds)] = 0
    # heatmap = gaussian_filter(heatmap, sigma=1.2)

        extent = [0, scen.sx, 0, scen.sy]
        # plt gt
        print(average_speeds.sum())
        scen.plot_value_field(axs, False, False)
        axs.imshow(average_speeds.T, extent=extent, origin="lower",
                   cmap="viridis", vmin=0, vmax=1.5)
        if export is not None:
            plt.savefig(export + ".svg", dpi=400)
        # plt.show()


    def plot_heatmap(self, axs, export=None, plot_example=False, t_max=-1,
                         binary=False, rw=False):
        x = []
        y = []
        t_vals = []
        source_x = []
        source_y = []
        map_name = self.experiments[0].field_name
        example_x = None
        example_y = None
        example_sources_x = None
        example_sources_y = None
        binsx = self.experiments[0].nx *2
        binsy = self.experiments[0].ny *2
        scen = ScenarioLoader(
            "../../resources_build/configs/maps/" + map_name + ".toml")
        scen.load_scenario()
        if binary:
            visited = np.zeros((binsx, binsy))
        for i, v in enumerate(self.experiments):
            t_vals = v.metrics["t"]
            if t_max == -1:
                t_max = t_vals[-1]
            t_id = np.argmin(np.abs(t_vals - t_max))
            x_v, y_v = v.get_states(rw=rw)
            x_v = x_v[:, :t_id]
            y_v = y_v[:, :t_id]
            print(t_id)
            print(t_vals[t_id])
            if plot_example and i == 7:
                example_x = x_v
                example_y = y_v
                example_sources_x, example_sources_y = v.get_sources("All")
            if binary:
                heatmap, xedges, yedges = np.histogram2d(x_v.flatten(),
                                                         y_v.flatten(),
                                                         bins=[binsx,
                                                               binsy],
                                                         range=[[0, scen.sx],
                                                                [0, scen.sy]])
                # add one to visited if that place has been visited at least once
                visited[np.where(heatmap > 0)] += 1

            x = x + x_v.flatten().tolist()
            y = y + y_v.flatten().tolist()

        x = np.array(x)
        y = np.array(y)
        if binary:
            heatmap = visited
        else:
            heatmap, xedges, yedges = np.histogram2d(x, y,
                                                     bins=[binsx,
                                                           binsy],
                                                     range=[[0, scen.sx],
                                                            [0, scen.sy]])
        # heatmap = gaussian_filter(heatmap, sigma=1.2)

        extent = [0, scen.sx, 0, scen.sy]
        colormap = cm.get_cmap("viridis")
        cmap = np.ones((colormap.N, 4))
        cmap[:, 1:3] = 0
        cmap[:, -1] = np.linspace(0.0, 1, colormap.N)
        cmap = ListedColormap(cmap)

        # plt gt
        scen.plot_value_field(axs, False, False)

        # imshow with colorbar
        plt.colorbar(axs.imshow(heatmap.T, extent=extent, origin="lower",
                   cmap=cmocean.cm.dense, interpolation="none"))

        if plot_example:
            axs.plot(example_x, example_y, c="white", linewidth=2)
            colors = cm.plasma(np.linspace(0, 1, len(example_x)))
            axs.scatter(example_x, example_y, c=colors)
            marker_style = dict(linestyle=':', marker='s',
                                markersize=100, markerfacecoloralt='tab:red')
            # axs.plot(example_sources_x, example_sources_y, c="white", linewidth=2)
            axs.scatter(example_sources_x, example_sources_y, c=colors,
                        marker="x", s=40)
            # last_id = [-1, -1]
            # id_counter = 0
            # label_id = 0
            # for i in range(len(example_x)):
            #     if last_id == [example_sources_x[i], example_sources_y[i]]:
            #         pass
            #     else:
            #         axs.annotate(label_id, (example_sources_x[i], example_sources_y[i]))
            #         label_id += 1
            #     last_id = [example_sources_x[i], example_sources_y[i]]

        # axs.scatter(sources[:, 0], sources[:, 1], c=cmap(norm(t)), s=100, marker="+")
        if export is not None:
            plt.savefig(export + ".svg", dpi=400)
        # plt.show()

    def get_metrics(self):
        self.metrics_keys = self.experiments[0].metrics_keys
        num_steps = self.experiments[0].num_steps
        for a in self.metrics_keys:
            self.metrics[a] = []
            for v in self.experiments:
                self.metrics[a].append(v.metrics[a])
                if len(v.metrics[a]) != num_steps and a != "source_threshold":
                    self.varying_lengths = True
                self.t_.append(v.metrics["t"])
        for a in self.metrics_keys:
            if not self.varying_lengths:
                met_arr = np.stack(self.metrics[a])
                self.mean_metrics[a] = met_arr.mean(axis=0)
                self.median_metrics[a] = np.median(met_arr, axis=0)
                self.std_metrics[a] = met_arr.std(axis=0)
                # self.std_metrics[a] = scipy.stats.iqr(met_arr, axis=0,interpolation = 'midpoint')
                self.q1_metrics[a] = np.quantile(met_arr, 0.25, axis=0)
                self.q2_metrics[a] = np.quantile(met_arr, 0.75, axis=0)
                self.t = self.experiments[0].metrics["t"]

    def plot_avg(self, axs, metric_keys, label="", first_time=None):
        t0 = 0
        if first_time is not None:
            # find entry in t that is closest to first_time
            t0 = np.abs(self.experiments[0].metrics["t"] - first_time).argmin()
            print(t0)

        for i, k in enumerate(metric_keys):
            if not self.varying_lengths:
                a = axs if len(metric_keys) == 1 else axs[i]
                a.plot(self.experiments[0].metrics["t"][t0::6],
                       self.median_metrics[k][t0::6],
                       label=label + self.name + (
                           "" if self.name == "" else "_") + k)
                # a.fill_between(self.experiments[0].metrics["t"][t0:],
                #                self.mean_metrics[k][t0:] - self.std_metrics[
                #                                                  k][
                #                                              t0:],
                #                self.mean_metrics[k][t0:] + self.std_metrics[
                #                                                  k][
                #                                              t0:], alpha=0.2)
                a.fill_between(self.experiments[0].metrics["t"][t0::6],
                               self.q1_metrics[k][t0::6],
                               self.q2_metrics[k][t0::6], alpha=0.2)

                a.set_xlabel("Time [s]")
                a.set_ylabel(k)
                print(self.median_metrics[k].mean())
            else:
                for j, v in enumerate(self.metrics[k]):
                    a = axs if len(metric_keys) == 1 else axs[i]
                    a.plot(self.t_[j], v, label=label + self.name + (
                        "" if self.name == "" else "_") + k)
                    a.plot(self.t_[j][-1], v[-1], "ok")

    def plot_source(self, axs):
        sources = []
        for a in self.experiments:
            x, y = a.get_sources(-1)
            sources.append([x, y])
        sources = np.stack(sources)
        # count how many times each source location is in the list
        unique, counts = np.unique(sources, axis=0, return_counts=True)
        # plot based on number of counter
        counts = counts / np.max(counts)
        axs.scatter(unique[:, 0], unique[:, 1], s=counts * 100, c="red",
                    alpha=0.7, marker="o")


class LogHandler(object):
    def __init__(self, path):
        self.path = path
        dirs = os.walk(self.path)
        a = next(dirs)
        self.is_single = False
        self.varying_params = False
        self.metrics_keys = None
        if len(a[2]) > 0:
            self.is_single = True
        else:
            b = next(dirs)
            if len(b[2]) > 0:
                self.varying_params = False
            else:
                self.varying_params = True
        self.experiments = None
        self.loadScenarioData()

    def plot_timestep(self, timestep, export=None):
        if self.is_single:
            self.experiments.plot_timestep(timestep, export)
        else:
            print("noo")

    def plot_final_histograms(self, export=None):
        # keys = ["CaseA", "CaseA1", "CaseA2", "CaseB", "CaseC"]
        keys = None
        if keys is None:
            keys = self.experiments.keys()
        colors = {
            "CaseA": mumred,
            "CaseB": mumgreen,
            "CaseC": mumblue,
            "CaseA1": mumred_alt,
            "CaseA2": mumred_alt
        }
        colors = None
        if not self.varying_params:
            print("Error: Cannot plot final histograms for single experiment")
            return
        errors = {}
        convergence_times = {}
        for k in keys:
            error, convergence_time = self.experiments[k].get_final_error()
            errors[self.experiments[k].name] = error
            convergence_times[self.experiments[k].name] = convergence_time
        fig, axs = plt.subplots()
        # axs.bar(errors.keys(), errors.values(), width=0.75, yerr=error_stds.values(), color=[colors[k] for k in errors.keys()], capsize=10)
        box_param = dict(patch_artist=True,
                         flierprops=dict(marker='.', markeredgecolor='black',
                                         fillstyle=None),
                         medianprops=dict(color='black'))
        boxplot = axs.boxplot(errors.values(),
                              positions=range(len(errors.keys())),
                              labels=errors.keys(), **box_param)
        if colors is not None:
            for i, name in enumerate(keys):
                boxplot['boxes'][i].set_facecolor(colors[name])
        axs.set_ylabel("Final Error")
        axs.grid()
        if not export is None:
            tikzplotlib.save(export + "_final_error.tex")
        plt.show()
        fig, axs = plt.subplots()

        boxplot = axs.boxplot(convergence_times.values(),
                              positions=range(len(convergence_times.keys())),
                              labels=convergence_times.keys(), **box_param)
        if colors is not None:
            for i, name in enumerate(keys):
                boxplot['boxes'][i].set_facecolor(colors[name])
        axs.set_ylabel("Convergence Time")
        axs.grid()
        if not export is None:
            tikzplotlib.save(export + "_convergence_time.tex")
        plt.show()

    '''
    Loads the scenario ground truth and simulation data
    '''

    def loadScenarioData(self):
        print("Is Single: ", self.is_single)
        print("Varying P: ", self.varying_params)
        if self.is_single:
            self.experiments = SingleExperiment(self.path)
            self.metrics_keys = self.experiments.metrics_keys
        elif not self.varying_params:
            self.experiments = RepeatedExperiment(self.path, "")
            self.metrics_keys = self.experiments.metrics_keys
        else:
            self.param_keys = []
            self.experiments = {}

            dirs = os.walk(self.path)
            a = next(dirs)
            for k in a[1]:
                self.param_keys.append(k)
                self.experiments[k] = RepeatedExperiment(
                    os.path.join(self.path, k), k)
            self.metrics_keys = self.experiments[
                self.param_keys[0]].metrics_keys

    def metric_overview(self, metric=None, export=False, ax=None, label=None,
                        first_time=None):
        metric_keys = []
        self_plot = False
        if metric is None:
            metric_keys = self.metrics_keys[1:]
        else:
            metric_keys.append(metric)
        if ax is None:
            self_plot = True
            fig, axs = plt.subplots(len(metric_keys), 1)
        else:
            axs = ax
            if len(metric_keys) == 1:
                if isinstance(axs, np.ndarray):
                    raise ValueError("Number of metrics and axes do not match")
            elif not axs.shape[0] == len(metric_keys):
                raise ValueError("Number of metrics and axes do not match")
        if self.is_single:
            for i, k in enumerate(metric_keys):
                axs[i].plot(self.experiments.metrics['t'],
                            self.experiments.metrics[k], label=label + k)
        elif not self.varying_params:
            self.experiments.plot_avg(axs, metric_keys, label,
                                      first_time=first_time)
        else:
            for i, a in enumerate(self.param_keys):
                self.experiments[a].plot_avg(axs, metric_keys, label,
                                             first_time=first_time)
        if len(metric_keys) == 1:
            pass
            # axs.legend()
        else:
            for a in range(axs.shape[0]):
                axs[a].legend()
        if export:
            tikzplotlib.save("metrics.tex")
        if self_plot:
            plt.show()

    def heatmap(self, export=None, plot_example=False):
        # plot heat map of all states
        fig, ax = plt.subplots()
        if self.is_single or not self.varying_params:
            self.experiments.plot_heatmap(ax, export, plot_example)
        else:
            for a in self.experiments:
                self.experiments[a].plot_heatmap(ax, export, plot_example)

    def plot_sources(self, axs=None, plot_map=False, export=None):
        gt_name = None
        if self.is_single or not self.varying_params:
            gt_name = self.experiments.field_name
        else:
            gt_name = self.experiments[self.param_keys[0]].field_name
        self_plot = False
        if axs is None:
            fig, axs = plt.subplots(1, 1)
            self_plot = True

        scen = ScenarioLoader(
            "../../resources_build/configs/maps/" + gt_name + ".toml")
        scen.load_scenario()
        scen.plot_value_field(axs, plot_map)
        scen.plot_quiver(axs)
        if self.is_single or not self.varying_params:
            self.experiments.plot_source(axs)
        else:
            for a in self.experiments:
                self.experiments[a].plot_source(axs)
        if export is not None:
            plt.savefig(export + ".svg", dpi=300)
        if self_plot:
            plt.show()

def get_first_convergence(data, threshold, tmax):
    # get's the first time step the data is consistently below the threshold for conv_length steps
    conv_length = 25
    # id_ = np.nan
    max_id = np.min([len(data), tmax])
    id_ = max_id
    for i in range(max_id):
        if np.all(data[i:i + conv_length] < threshold):
            id_ = i
            break
    return id_
def loadMetricDataFrame(path, subdirs=None, t=-1, cases=None, paramvals=None,
                        get_abs=False):
    if subdirs is None:
        subdirs = os.listdir(path)
    data = {}
    data["Experiment"] = []
    data["Case"] = []
    data["ParamVal"] = []
    # treshold DICT todo make this automatic
    thresholds = {}
    thresholds["CaseA"] = 1.0
    thresholds["CaseB"] = 5.0
    thresholds["CaseC"] = 2.0
    thresholds["CaseA1"] = 1.5
    thresholds["CaseA2"] = 2.0

    for directory in subdirs:
        exp_dir = os.path.join(path, directory)
        case_dirs = os.listdir(exp_dir) if cases is None else cases
        for case in case_dirs:
            case_dir = os.path.join(exp_dir, case)
            param_dirs = os.listdir(
                case_dir) if paramvals is None else paramvals
            for param in param_dirs:
                param_dir = os.path.join(case_dir, param)
                for repeat in os.listdir(param_dir):
                    repeat_dir = os.path.join(param_dir, repeat)
                    exp = SingleExperiment(repeat_dir)
                    data["Experiment"].append(directory)
                    data["Case"].append(case)
                    data["ParamVal"].append(param)
                    metrics = exp.metrics
                    if t == -1:
                        t_id = -1
                    else:
                        t_id = np.argmin(np.abs(np.array(metrics["t"]) - t))
                    for k in metrics.keys():
                        if not k in data.keys():
                            data[k] = []
                        data[k].append(metrics[k][t_id])
                    # get when source dist falls under 1.0
                    k = "SourceDistMetric"
                    k_threshold = "SourceDistThreshold"
                    if not k_threshold in data.keys():
                        data[k_threshold] = []
                    # get first time when source dist falls under 1.0

                    id_ = get_first_convergence(metrics[k], thresholds[case], t_id)
                    if np.isnan(id_):
                        data[k_threshold].append(np.nan)
                    elif id_ > t_id:
                        data[k_threshold].append(metrics["t"][t_id])
                    else:
                        data[k_threshold].append(metrics["t"][id_])
                    if get_abs:
                        # now add the abs error in the field
                        k = "abs_error"
                        if not k in data.keys():
                            data[k] = []
                        data[k].append(exp.getErrorMap(t_id))
    return pd.DataFrame(data)


if __name__ == "__main__":
    # lh = LogHandler("../../Experiments/Paper/Random2/RoomFieldTest3/")
    # lh2 = LogHandler("../../Experiments/Paper/Random2/RoomFieldTest4/")
    # case = "CaseB"
    # num_nodes = "40"
    fig, axs = plt.subplots(1, 1)
    s = SingleExperiment(
        "../../Experiments/FinalEval/FlowTHESIS/CaseB/50/0/")
    s.read_data()
    s.plot_mean_error(axs)
    plt.show()
    # s2 = SingleExperiment(
    #     "../../Experiments/FinalEval/NoFlow2/" + case + "/" + num_nodes + "/0/")
    # lh = LogHandler(
    #     "../../Experiments/FinalEval/NoFlow/" + case + "/" + num_nodes + "/")
    # lh2 = LogHandler(
    #     "../../Experiments/FinalEval/Flow/" + case + "/" + num_nodes + "/")
    # fig, axs = plt.subplots(4, 1)
    # lh.metric_overview(ax=axs, label="NoFlow")
    # # fig, axs = plt.subplots(4, 1)
    # lh2.metric_overview(ax=axs, label="Flow")
    # plt.show()
    # # s = SingleExperiment("../../Experiments/Paper/Concept/0/")
    # # s2 = SingleExperiment("../../Experiments/Paper/Concept/1/")
    # fig, axs = plt.subplots(3, 5)
    # s.plot_evaluation(axs[:, 0:2], t=-1)
    # s2.plot_evaluation(axs[:, 2:4], t=-1)
    # gt = ScenarioLoader("../../resources_build/configs/maps/CaseA1.toml")
    # gt.load_scenario()
    # gt.plot_value_field(axs[0, 4], plot_values=True)
    #
    # plt.show()
    # plt.savefig("concept.svg", dpi=300)
    # plt.show()
    #
    # ts = 20
    # lh.plot_timestep(ts, export=ts)
    # lh_bfs = LogHandler("../../Experiments/Paper/BFS2/")
    # lh_bfs.plot_final_histograms()
    # lh = LogHandler("../../Experiments/Paper/Random4/CaseA1/")
    # lh2 = LogHandler("../../Experiments/Paper/Random4/CaseA/")
    # lh3 = LogHandler("../../Experiments/Paper/Random4/CaseB/")
    # lh4 = LogHandler("../../Experiments/Paper/Random4/CaseC/")
    # lh5 = LogHandler("../../Experiments/Paper/Random4/CaseA2/")
    # lh = LogHandler("../../Experiments/FinalEval/Random/CaseA1/50/")
    # lh2 = LogHandler("../../Experiments/FinalEval/Random/CaseA/50/")
    # lh3 = LogHandler("../../Experiments/FinalEval/Random/CaseB/50/")
    # lh4 = LogHandler("../../Experiments/FinalEval/Random/CaseC/50/")
    # lh5 = LogHandler("../../Experiments/FinalEval/Random/CaseA2/50/")
    # lh2 = LogHandler("../../Experiments/Paper/BFS2/RoomFieldTest4")
    # lh2 = LogHandler("../../Experiments/avg100/RoomField2")
    # lh = LogHandler("../../Experiments/BFS/RoomFieldC/")
    # lh2 = LogHandler("../../Experiments/BFS/RoomField2/")
    # # lh = LogHandler("../../Experiments/avg100/RoomField2/")
    # fig, axs = plt.subplots(4, 1)
    # lh.metric_overview(ax=axs, label="Random")
    # # fig, axs = plt.subplots(4, 1)
    # lh2.metric_overview(ax=axs, label="noFlow")
    # plt.show()
    # plot_boxplot_metric("RMSEMetric", -1, "../../Experiments/ClosedLoop/NumberOfAgents"
    #                                       "/Exploration/")
    fig, axs = plt.subplots(1, 1)
    lh1 = LogHandler(
        "../../Experiments/ClosedLoop/NumberOfAgents/SourceLoc7/CaseA/1/")
    lh2 = LogHandler(
        "../../Experiments/ClosedLoop/NumberOfAgents/SourceLoc7/CaseA/2/")
    lh3 = LogHandler(
        "../../Experiments/ClosedLoop/NumberOfAgents/SourceLoc7/CaseA/3/")
    lh4 = LogHandler(
        "../../Experiments/ClosedLoop/NumberOfAgents/SourceLoc7/CaseA/4/")
    # lh1 = LogHandler(
    #     "../../Experiments/ClosedLoop/NumberOfAgents/ExplorationNewParams/CaseA/4/")
    # lh2 = LogHandler(
    #     "../../Experiments/ClosedLoop/NumberOfAgents/ExplorationBallistic2/CaseA/4/")
    # lh3 = LogHandler(
    #     "../../Experiments/ClosedLoop/NumberOfAgents/ExplorationBrownian2/CaseA/4/")
    # # lh3 = LogHandler(
    # #     "../../Experiments/FinalEval/FlowFinal/CaseA/50/")
    # # lh4 = LogHandler(
    # #     "../../Experiments/FinalEval/FLowNew/CaseA/50/")
    # # lh1_rw = LogHandler("../../Experiments/ClosedLoop/NumberOfAgents/ExplorationBallistic/CaseA2/1/")
    # # lh2_rw = LogHandler("../../Experiments/ClosedLoop/NumberOfAgents/ExplorationBallistic/CaseA2/2/")
    # # lh3_rw = LogHandler("../../Experiments/ClosedLoop/NumberOfAgents/ExplorationBallistic/CaseA2/3/")
    # # lh4_rw = LogHandler("../../Experiments/ClosedLoop/NumberOfAgents/ExplorationBallistic/CaseA2/4/")
    # # lh1.metric_overview("FieldLogLikelihoodMetric", ax=axs, label="30")
    lh1.metric_overview("RMSEMetric", ax=axs, label="1")
    lh2.metric_overview("RMSEMetric", ax=axs, label="2")
    lh3.metric_overview("RMSEMetric", ax=axs, label="3")
    lh4.metric_overview("RMSEMetric", ax=axs, label="4")
    axs.grid()
    axs.legend()
    tikzplotlib.save("rmse_sourcedist_A.tex")
    # axs.grid()
    # axs.set_xlabel("Number of Measurements")
    # # lh4.metric_overview("RMSEMetric", ax=axs, label="4")
    # # tikzplotlib.save("sourceDist_caseA.tex")
    # # lh1_rw.metric_overview(ax=axs[:, 0], label="1_RW")
    # # lh2_rw.metric_overview(ax=axs[:, 0], label="2_RW")
    # # lh3_rw.metric_overview(ax=axs[:, 0], label="3_RW")
    # # lh4_rw.metric_overview(ax=axs[:, 0], label="4_RW")
    # tikzplotlib.save("caseA_ll_flowvsnoflow.tex")
    # plt.show()
    # single_flow = SingleExperiment(
    #     "../../Experiments/FinalEval/FlowFinal/CaseB/50/1/")
    # single_flow.read_data()
    # fig, axs = plt.subplots(1,1)
    #
    # single_flow.plot_mean_error(axs)
    # plt.savefig("errorCaseB.svg")
    # plt.show()
    # fig, axs = plt.subplots(1, 1)
    #
    # single_flow.plot_estimates(axs, t=-1)
    # scen = ScenarioLoader("../../resources_build/configs/maps/CaseA1.toml")
    # scen.load_scenario()
    # scen.plot_value_field(axs, plot_values=False)
    #
    plt.show()
    # lh1 = LogHandler("../../Experiments/ClosedLoop/ExampleRuns/")
    # lh1.heatmap()

    # s = LogHandler("../../Experiments/FinalEval/FlowTHESIS/CaseB/50/")
    # fig, axs = plt.subplots(1, 1)
    # s.plot_sources(axs=axs, plot_map=True)
    # plt.savefig("sources_caseB_50.svg")
    # plt.show()
    # lh.heatmap("heatmap", plot_example=True)
