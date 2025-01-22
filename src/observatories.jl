using AstroLib

AstroLib.observatories["icecube"] = AstroLib.Observatory(
    "IceCube",
    -90 * units.degree,
    0 * units.degree,
    -1947 * units.m,
    0.0
)

AstroLib.observatories["pone"] = AstroLib.Observatory(
    "P-ONE",
    48.4284 * units.degree,
    -123.3656 * units.degree,
    -1947 * units.m,
    -8
)

AstroLib.observatories["arca"] = AstroLib.Observatory(
    "KM3NeT-ARCA",
    ten(36, 16) * units.degree,
    ten(16, 6) * units.degree,
    -1947 * units.m,
    1
)

AstroLib.observatories["orca"] = AstroLib.Observatory(
    "KM3NeT-ORCA",
    ten(42, 48) * units.degree,
    ten(6, 2) * units.degree,
    -1947 * units.m,
    1
)

AstroLib.observatories["gvd"] = AstroLib.Observatory(
    "Baikal-GVD",
    53.5587*units.degree,
    108.1650*units.degree, 
    -800units.m,
    -7
)

AstroLib.observatories["tambo"] = AstroLib.Observatory(
    "TAMBO",
    larger_valley_coord.latitude,
    larger_valley_coord.longitude,
    2units.m,
    -7
)

AstroLib.observatories["trident"] = AstroLib.Observatory(
    "TRIDENT",
    17.4 * units.degree,
    114.0 * units.degree,
    -3units.km,
    -7
)

AstroLib.observatories["eq"] = AstroLib.Observatory(
    "EQ",
    0.0 * units.degree,
    0.0 * units.degree,
    2units.m,
    -7
)
