import cdspyreadme

tablemaker = cdspyreadme.CDSTablesMaker()

table1 = tablemaker.addTable("./data/toi4438_spectroscopy.csv", name="tablec1.dat",
                              description="TOI-4438 spectroscopic time series")
table2 = tablemaker.addTable("./data/toi442_spectroscopy.csv", name="tablec2.dat",
                              description="TOI-442 spectroscopic time series")
table3 = tablemaker.addTable("./data/toi4438_polarimetry.csv", name="tablec3.dat",
                              description="TOI-4438 polarimetric time series")
table4 = tablemaker.addTable("./data/toi442_polarimetry.csv", name="tablec4.dat",
                              description="TOI-442 polarimetric time series")

spec_columns = {
    "rjd":       ("d",     "Reduced Julian Date (JD-2400000)",   "F16.11"),
    "rv":        ("m/s",   "Radial velocity",                                     "F9.1"),
    "rv_err":    ("m/s",   "Radial velocity uncertainty",                         "F5.1"),
    "dlW":       ("m2/s2", "Differential line width activity indicator",          "F10.0"),
    "dlW_err":   ("m2/s2", "Differential line width uncertainty",                 "F7.0"),
    "dtemp":     ("K",     "Line-depth differential temperature activity indicator", "F6.1"),
    "dtemp_err": ("K",     "Differential temperature indicator uncertainty",      "F5.1"),
}
pol_columns = {
    "rjd": ("d", "Reduced Julian Date (JD-2400000)", "F12.6"),
    "bl":  ("G", "Longitudinal magnetic field",                        "F6.1"),
    "sbl": ("G", "Longitudinal magnetic field uncertainty",            "F5.1"),
}

for table in (table1, table2):
    for name, (unit, description, fmt) in spec_columns.items():
        col = table.get_column(name)
        col.unit = unit
        col.description = description
        col.set_format(fmt)

for table in (table3, table4):
    for name, (unit, description, fmt) in pol_columns.items():
        col = table.get_column(name)
        col.unit = unit
        col.description = description
        col.set_format(fmt)

# write table in CDS-ASCII aligned format (required)
tablemaker.writeCDSTables()

# Customize ReadMe output
tablemaker.title = "Revisiting TOI-4438 and TOI-442 planetary systems with new observations from SPIRou and TESS"
tablemaker.author = "J. Serrano Bell"
tablemaker.authors = "J. Serrano Bell, G. Hébrard, E. Martioli, R. F. Díaz, L. de Almeida, D. Lorenzo-Oliveira, A. Salmi, C. Dorn, M. Valatsou, A. Carmona, M. Ould-Elhkim, L. Arnold, É. Artigau, I. Boisse, X. Bonfils, C. Cadieux, Z. Chakir, N. J. Cook, X. Delfosse, J.-F. Donati, R. Doyon, N. Heidari, J. M. Jenkins, F. Kiefer, S. Lafrance, A. L'Heureux, C. Moutou, J. Morneau, X. Vandelac"
tablemaker.date = 2026
tablemaker.abstract = "We present a comprehensive re-analysis of two star-planet systems: TOI-4438, an M3.5V star hosting a mini-Neptune in a 7.4-day orbit, and TOI-442, an M1V star with a 4-day period planet located within the hot Neptune desert. Both systems were originally identified as transiting planet candidates by the Transiting Exoplanet Survey Satellite (TESS) and subsequently validated through the radial velocity (RV) method. Our work incorporates new TESS transit data and high-resolution spectroscopy from the SPIRou near-infrared (nIR) spectropolarimeter. We detect a persistent and relatively strong Zeeman signature in TOI-442, while TOI-4438 exhibits weaker and intermittent magnetic activity, and we infer the stellar rotation periods of both stars from the variability of the longitudinal magnetic field. We jointly fit photometry and RV models for each system. For TOI-4438 b we combine archival CARMENES data with 81 SPIRou observations and five TESS sectors. This yields a refined planetary mass of Mp = 4.11(+0.40)(-0.38)Mearth and a radius of Rp = 2.40(+0.09)(-0.10)Rearth, consistent with the previous estimate within 1.1 sigma while improving by 53% the precision on the mass and 22% on the radius. For TOI-442 b, we add 29 SPIRou RV measurements to an extensive archival dataset, significantly extending the temporal baseline. Incorporating this with a new TESS sector, we tighten the constraints on the planetary mass to Mp = 28.38(+0.77)(-0.73)Mearth and radius to Rp = 4.25(+0.10)(-0.08)Rearth, which agrees to the previous values within 1.5 sigma and improves the precisions by 46% and 67% respectively. We find no clear signs of additional planets in the currently available RV data, although we detect a single-transit event in the TOI-4438 light curve. We compare various RV models and find that those accounting for stellar variability-induced signals yield improved constraints on the planetary parameters."

tablemaker.license = "CC-BY-4.0"
tablemaker.more_description = ("Radial velocities and spectral activity indicators (differential line "
    "width, line-depth temperature proxy) were extracted from SPIRou spectra using the line-by-line "
    "(LBL) method. Longitudinal magnetic field measurements were derived from the circularly polarized "
    "(Stokes V) spectra using Least-Squares Deconvolution.")  # TODO: confirm pipeline details/references
tablemaker.keywords = "planets and satellites: detection, techniques: photometric, techniques: radial velocities, stars: planetary systems"
# TODO: add once available, e.g. tablemaker.bibcode = "2026A&A...676A...1S"
# tablemaker.putRef("II/246", "2mass catalogue")
# tablemaker.putRef("http://...", "external link")
tablemaker.acknowledgements = "Juan Serrano Bell, jserrano(at)unsam.edu.ar"


def _object_line(name, ra=None, de=None):
    """ra/de as copied straight from SIMBAD: 'hh mm ss.ss' / '+dd mm ss.ss'."""
    if ra is None or de is None:
        return f"    {'TODO RA':<11}{'TODO DE':<10}{name}"
    return f"    {ra}  {de}    {name}"


# "Objects:", "Acknowledgements:" and a null "See also:" are static/broken in cdspyreadme's
# packaged ReadMe.template (no $objects var, and printRef() returns None with no putRef()
# calls, which leaks the literal "None" into the file). ./ReadMe.template is a copy of the
# packaged template with $objects/$acknowledgements added and $seealso guarded.
objects = "\n".join([
    _object_line("TOI-4438", "18 01 16.10", "+35 35 50.49"), 
    _object_line("TOI-442", "04 16 45.60", "-12 05 02.45"),  
])
_seealso = tablemaker.printRef(outBuffer=True) or ""
tablemaker.setReadmeTemplate(
    "./ReadMe.template",
    templateValue={"objects": objects, "acknowledgements": tablemaker.acknowledgements, "seealso": _seealso},
)

# Save ReadMe (in current directory)
tablemaker.makeReadMe()

# print ReadMe
import sys
tablemaker.makeReadMe(out=sys.stdout)