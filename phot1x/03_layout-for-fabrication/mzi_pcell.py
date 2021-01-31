from siepicfab import all as pdk
from ipkiss3 import all as i3
from ipkiss.technology import get_technology
from circuit.all import CircuitCell
from circuit.connector_functions import manhattan
from functools import partial
import pylab as plt
import numpy as np

TECH = get_technology()


class MZI(CircuitCell):

    control_points = i3.ListProperty(doc="Position of the control point for the longer arm of the MZI")
    bend_radius = i3.PositiveNumberProperty(default=5.0, doc="Bend radius of the waveguides")

    fgc = i3.ChildCellProperty(locked=True, doc="PCell for the fiber grating coupler")
    splitter = i3.ChildCellProperty(locked=True, doc="PCell for the Y-Branch")
    dir_coupler = i3.ChildCellProperty(locked=True, doc="PCell for the directional coupler")

    def _default_control_points(self):
        return [(100.0, 220.0)]

    def _default_fgc(self):
        return pdk.EbeamGCTE1550()

    def _default_splitter(self):
        return pdk.EbeamY1550()

    def _default_dir_coupler(self):
        return pdk.EbeamBDCTE1550()

    def _default_child_cells(self):
        child_cells = {
            "fgc_1": self.fgc,
            "fgc_2": self.fgc,
            "fgc_3": self.fgc,
            "yb": self.splitter,
            "dc": self.dir_coupler,
        }
        return child_cells

    def _default_joins(self):
        joins = [("fgc_3:opt1", "yb:opt1")]
        return joins

    def _default_connectors(self):
        br = self.bend_radius
        connectors = [
            ("yb:opt3", "dc:opt2", manhattan, {"bend_radius": br}),
            ("yb:opt2", "dc:opt1", partial(manhattan, control_points=self.control_points), {"bend_radius": br}),
            ("dc:opt4", "fgc_2:opt1", manhattan, {"bend_radius": br}),
            ("dc:opt3", "fgc_1:opt1", manhattan, {"bend_radius": br}),
        ]
        return connectors

    def _default_place_specs(self):
        place_specs = [
            i3.Place("fgc_1:opt1", (0, 0)),
            i3.PlaceRelative("fgc_2:opt1", "fgc_1:opt1", (0.0, 120.0)),
            i3.PlaceRelative("fgc_3:opt1", "fgc_2:opt1", (0.0, 120.0)),
            i3.PlaceRelative("dc:opt1", "yb:opt2", (20.0, -40.0), angle=90),
        ]
        return place_specs

    def _default_external_port_names(self):
        epn = {
            "fgc_3:fib1": "in",
            "fgc_2:fib1": "out1",
            "fgc_1:fib1": "out2",
        }
        return epn


if __name__ == "__main__":

    # Layout
    mzi = mzi = MZI(
        control_points=[(100.0, 220.0)],
        bend_radius=5.0,
    )
    mzi_layout = mzi.Layout()
    mzi_layout.visualize(annotate=True)
    mzi_layout.visualize_2d()
    mzi_layout.write_gdsii("mzi.gds")

    # Circuit model
    my_circuit_cm = mzi.CircuitModel()
    wavelengths = np.linspace(1.5, 1.6, 10001)
    S_total = my_circuit_cm.get_smatrix(wavelengths=wavelengths)

    # Plotting
    plt.plot(wavelengths, 20 * np.log10(np.abs(S_total['out1:0', 'in:0'])), '-', linewidth=2.2, label="TE-out1")
    plt.plot(wavelengths, 20 * np.log10(np.abs(S_total['out2:0', 'in:0'])), '-', linewidth=2.2, label="TE-out2")
    plt.xlabel('Wavelength [um]', fontsize=16)
    plt.ylabel('Transmission [dB]', fontsize=16)
    plt.legend(fontsize=14, loc=4)
    plt.show()
