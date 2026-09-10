import json
import math
from dataclasses import replace

import numpy as np
import pytest

from sminus2_point_particle.config import EvolutionConfig, OutputRequest
from sminus2_point_particle.errors import ContractError, OpenBoundaryError
from sminus2_point_particle.evolution import continuum_rhs
from sminus2_point_particle.extraction import project_mode, projection_row, scri_slice, spin_weighted_spherical_harmonic
from sminus2_point_particle.field import angular_operator, coefficients, p_from_pi
from sminus2_point_particle.geometry import bl_timelike_norm, boyer_lindquist_radius, compactified_radius, isco_quantities, radial_state
from sminus2_point_particle.initial_data import FluxRecord, generate_ori_thorne_initial_data, initial_data_from_manifest
from sminus2_point_particle.operators import BoundaryRowExtension, interior_derivative, open_uniform_nodes
from sminus2_point_particle.source import evaluate_source_blocks
from sminus2_point_particle.trajectory import stage_sample


def test_config_round_trip_normalization_and_stable_key():
    config = EvolutionConfig(62.0, 0.7, 2, 64, 32, 100.0, flux_model_version="flux-a")
    restored = EvolutionConfig.from_json(config.to_json())
    assert restored == config
    assert restored.M == restored.L == 1.0
    assert restored.q_mass == 1e-5
    assert config.initial_data_key("1", "abc") == restored.initial_data_key("1", "abc")
    assert json.loads(config.to_json())["mass_scale"] == 62.0
    assert json.loads(config.to_json())["q_mass"] == 1e-5
    output = OutputRequest((2, 4), save_scri_field=True)
    assert OutputRequest.from_json(output.to_json()) == output
    with pytest.raises(ContractError):
        EvolutionConfig(62.0, 0.7, 2, 64, 32, 100.0, transition_prescription="other")
    with pytest.raises(ContractError):
        EvolutionConfig(62.0, 0.7, 2, 64, 32, 100.0).initial_data_key("1", "abc")


def test_schwarzschild_isco_and_versioned_transition_manifest():
    r, energy, angular_momentum, omega = isco_quantities(0.0)
    assert r == pytest.approx(6.0)
    assert energy == pytest.approx(math.sqrt(8.0 / 9.0))
    assert angular_momentum == pytest.approx(2 * math.sqrt(3.0))
    assert omega == pytest.approx(1 / (6 * math.sqrt(6.0)))
    plunge = generate_ori_thorne_initial_data(FluxRecord(0.0, 6, 0.001, 0.01, "flux-test-v1", "unit fixture"))
    manifest = plunge.manifest
    assert plunge.r0 == pytest.approx(r - manifest.delta_r)
    with pytest.raises(ContractError):
        initial_data_from_manifest(1.0, 0.0, 1.0, replace(manifest, prescription="other"))
    with pytest.raises(ContractError):
        initial_data_from_manifest(1.0, 0.0, 1.0, replace(manifest, q_mass=1e-3))


def test_radius_inverse_and_velocity_contracts():
    r = np.array([2.5, 4.0, 9.0])
    assert np.allclose(boyer_lindquist_radius(compactified_radius(r, 1.3), 1.3), r)
    state = radial_state(5.0, 1.0, 0.5, 1.0, 0.95, 3.0)
    assert bl_timelike_norm(state, 1.0, 0.5) == pytest.approx(-1.0, abs=2e-14)
    # Independent coordinate-Jacobian expressions from DERIVATION.md Section III.
    de = 25.0 - 10.0 + 0.25
    hprime = (25.0 + 0.25) / de - 2.0 - 4.0 / 5.0
    assert state.u_T == pytest.approx(state.u_t_bl + hprime * state.u_r_bl)
    assert state.u_R == pytest.approx(-state.u_r_bl / 25.0)
    assert state.u_Phi == pytest.approx(state.u_phi_bl + 0.5 * state.u_r_bl / de)
    h = 1e-3
    plus = radial_state(5.0 + h, 1.0, 0.5, 1.0, 0.95, 3.0)
    minus = radial_state(5.0 - h, 1.0, 0.5, 1.0, 0.95, 3.0)
    for value, first, second, value_plus, value_minus in (
        (state.u_T, state.u_T_r, state.u_T_rr, plus.u_T, minus.u_T),
        (state.u_R, state.u_R_r, state.u_R_rr, plus.u_R, minus.u_R),
        (state.u_Phi, state.u_Phi_r, state.u_Phi_rr, plus.u_Phi, minus.u_Phi),
    ):
        assert (value_plus - value_minus) / (2 * h) == pytest.approx(first, abs=2e-7)
        assert (value_plus - 2 * value + value_minus) / h**2 == pytest.approx(second, abs=2e-7)


def test_stage_sample_trajectory_and_mass_linearity():
    s1 = stage_sample(0.0, 5.0, 0.3, 1.0, 0.5, 1.0, 0.95, 3.0, particle_mass=1.0)
    s3 = stage_sample(0.0, 5.0, 0.3, 1.0, 0.5, 1.0, 0.95, 3.0, particle_mass=3.0)
    assert s1.R_jet.first == pytest.approx(s1.u_R.value / s1.u_T.value)
    assert s1.Phi_jet.first == pytest.approx(s1.u_Phi.value / s1.u_T.value)
    assert s3.amplitude_nn.value == pytest.approx(3 * s1.amplitude_nn.value)
    assert s3.amplitude_mn.second == pytest.approx(3 * s1.amplitude_mn.second)
    assert s3.amplitude_mm.first == pytest.approx(3 * s1.amplitude_mm.first)
    with pytest.raises(ContractError):
        stage_sample(0.0, 5.0, 0.3, 1.0, 0.5, 1.0, 0.95, 3.0, particle_mass=-1.0)


def test_interior_operators_are_exact_on_polynomials_and_boundaries_stay_open():
    x = np.linspace(0.0, 1.0, 13)
    d1 = interior_derivative(x, 1, 5)
    d2 = interior_derivative(x, 2, 5)
    f = x**4 - 2 * x**3 + x
    assert np.allclose(d1.apply(f)[d1.valid_rows], (4 * x**3 - 6 * x**2 + 1)[d1.valid_rows], atol=2e-13)
    assert np.allclose(d2.apply(f)[d2.valid_rows], (12 * x**2 - 12 * x)[d2.valid_rows], atol=2e-12)
    assert np.isnan(d1.apply(f)[0]) and np.isnan(d1.apply(f)[-1])
    with pytest.raises(OpenBoundaryError):
        BoundaryRowExtension().require_production()


def test_field_coefficients_against_independent_reference_and_p_relation():
    R, y, M, a, L, m = 0.17, -0.23, 1.0, 0.6, 1.0, 2
    got = coefficients(R, y, M, a, L, m)
    refs = (
        8 * M * (2 * M - a * a * R / L**2) * (1 + 2 * M * R / L**2) - a * a * (1 - y * y),
        -2 * (L**2 - (8 * M * M - a * a) * R**2 / L**2 + 4 * a * a * M * R**3 / L**4),
        -(L**2 - 2 * M * R + a * a * R**2 / L**2) * R**2 / L**2,
        2j * a * m * (1 + 4 * M * R / L**2) + 2 * (2 * M * (2 - 3 * a * a * R**2 / L**4) - a * a * R / L**2 + 2j * a * y),
        2j * a * m * R**2 / L**2 + 2 * R * (1 + M * R / L**2 - 2 * a * a * R**2 / L**4),
        2j * a * m * R / L**2 - 2 * M * R / L**2 - 2 * a * a * R**2 / L**4,
    )
    for value, ref in zip((got.A, got.B, got.C, got.D, got.E_R, got.F), refs):
        assert value == pytest.approx(ref)
    psi, pi, q = 1.2 - 0.4j, -0.7 + 0.3j, 0.2 + 0.8j
    P = p_from_pi(psi, pi, q, got)
    rhs = continuum_rhs(psi, P, q, -0.11j, 0.9 + 0.1j, 0.4 - 0.2j, got)
    assert rhs.dpsi == pytest.approx(pi)
    expected_dP = (0.4 - 0.2j) - got.C * (-0.11j) + (0.9 + 0.1j) - got.E_R * q - got.F * psi
    assert rhs.dP == pytest.approx(expected_dP)


def test_angular_operator_independent_expression_and_axis_rejection():
    y = np.array([-0.4, 0.1, 0.7])
    psi = 1 + y + 2 * y**2
    got = angular_operator(psi, 1 + 4 * y, np.full_like(y, 4.0), y, 2)
    ref = (1 - y**2) * 4 - 2 * y * (1 + 4 * y) - (2 + 2 * y) ** 2 * psi / (1 - y**2) - 2 * psi
    assert np.allclose(got, ref)
    with pytest.raises(ContractError):
        angular_operator(1.0, 0.0, 0.0, 1.0, 2)


def test_four_source_blocks_are_separate_sum_consistent_and_mu_linear():
    rnodes = open_uniform_nodes(0.06, 0.46, 25)
    ynodes = open_uniform_nodes(-1.0, 1.0, 25)
    R, y = np.meshgrid(rnodes, ynodes, indexing="ij")
    dr = interior_derivative(rnodes, 1, 5)
    dy = interior_derivative(ynodes, 1, 5)
    derivative_r = lambda value: dr.apply(value, axis=0)
    derivative_y = lambda value: dy.apply(value, axis=1)
    s1 = stage_sample(0.0, 5.0, 0.2, 1.0, 0.5, 1.0, 0.95, 3.0, particle_mass=1.0)
    s2 = stage_sample(0.0, 5.0, 0.2, 1.0, 0.5, 1.0, 0.95, 3.0, particle_mass=2.0)
    b1 = evaluate_source_blocks(R, y, s1, 1.0, 0.5, 1.0, 2, 0.035, 0.12, derivative_r, derivative_y)
    b2 = evaluate_source_blocks(R, y, s2, 1.0, 0.5, 1.0, 2, 0.035, 0.12, derivative_r, derivative_y)
    core = np.s_[4:-4, 4:-4]
    manual_total = b1.block1 + b1.block2 + b1.block3 + b1.block4
    assert np.allclose(b1.total[core], manual_total[core], rtol=0, atol=0)
    assert np.all(np.isfinite(b1.total[core]))
    for first, second in zip((b1.block1, b1.block2, b1.block3, b1.block4), (b2.block1, b2.block2, b2.block3, b2.block4)):
        assert np.allclose(second[core], 2 * first[core], rtol=2e-12, atol=1e-15)


def test_four_source_block_formulas_with_independent_zero_spatial_derivative_reference():
    R, y = np.array([[0.2]]), np.array([[-0.15]])
    M, a, L, m = 1.0, 0.5, 1.0, 2
    sample = stage_sample(0.0, 5.0, 0.2, M, a, L, 0.95, 3.0)
    zero = lambda value: np.zeros_like(value)
    got = evaluate_source_blocks(R, y, sample, M, a, L, m, 0.04, 0.13, zero, zero)
    x = R - sample.R
    gaussian = np.exp(-x**2 / (2 * 0.04**2) - y**2 / (2 * 0.13**2)) / (2 * np.pi * 0.04 * 0.13)
    gamma1 = x * sample.R_jet.first / 0.04**2
    gamma2 = x**2 * sample.R_jet.first**2 / 0.04**4 + (x * sample.R_jet.second - sample.R_jet.first**2) / 0.04**2
    def dressed(amplitude):
        phase = np.exp(-1j * m * sample.Phi)
        q0 = phase * amplitude.value
        q1 = phase * (amplitude.first - 1j * m * sample.Phi_jet.first * amplitude.value)
        q2 = phase * (amplitude.second - 2j * m * sample.Phi_jet.first * amplitude.first - 1j * m * sample.Phi_jet.second * amplitude.value - m**2 * sample.Phi_jet.first**2 * amplitude.value)
        return q0 * gaussian, (q1 + q0 * gamma1) * gaussian, (q2 + 2 * q1 * gamma1 + q0 * gamma2) * gaussian
    fnn, fmn, fmm = dressed(sample.amplitude_nn), dressed(sample.amplitude_mn), dressed(sample.amplitude_mm)
    rho = R / (L**2 + 1j * a * R * y)
    rb = R / (L**2 - 1j * a * R * y)
    c = np.sqrt(1 - y**2)
    ell0, ellm1 = m / c, (m + y) / c
    jt, jr = -(2 + 4 * M * R / L**2), -R**2 / L**2
    pref = -16 * np.pi * R * (L**4 + a**2 * R**2 * y**2)
    u10, u11, u12 = [rho**-2 * rb**-1 * f for f in fnn]
    v10, v11 = -1j * a * c * u11 + ell0 * u10, -1j * a * c * u12 + ell0 * u11
    w10, w11 = rho**-4 * v10, rho**-4 * v11
    ref1 = pref * (-rho**8 * rb * (-1j * a * c * w11 + ellm1 * w10))
    u20, u21, u22 = [rho**-2 * rb**-2 * R**2 * f for f in fmn]
    w20, w21 = rho**-4 * rb**2 * (jt * u21), rho**-4 * rb**2 * (jt * u22)
    ref2 = pref * (-rho**8 * rb / (np.sqrt(2) * R**4) * (-1j * a * c * w21 + ellm1 * w20))
    u30, u31, u32 = [rho**-2 * rb * f for f in fmm]
    w30, w31 = rho**-4 * (jt * u31), rho**-4 * (jt * u32)
    ref3 = pref * (-rho**8 * rb / (2 * R**4) * (jt * w31))
    u40, u41, u42 = [rho**-2 * rb**-2 * f for f in fmn]
    v40, v41 = -1j * a * c * u41 + ellm1 * u40, -1j * a * c * u42 + ellm1 * u41
    w40, w41 = rho**-4 * rb**2 * R**2 * v40, rho**-4 * rb**2 * R**2 * v41
    ref4 = pref * (-rho**8 * rb / (np.sqrt(2) * R**4) * (jt * w41))
    for value, reference in zip((got.block1, got.block2, got.block3, got.block4), (ref1, ref2, ref3, ref4)):
        assert np.allclose(value, reference, rtol=2e-14, atol=1e-16)


def test_scri_projection_recovers_analytic_22_mode_and_normalization():
    y, weights = np.polynomial.legendre.leggauss(48)
    harmonic = spin_weighted_spherical_harmonic(2, 2, y)
    # Independent closed form in y=-cos(theta): _{-2}Y_22=sqrt(5/(64*pi))*(1-y)^2.
    explicit = math.sqrt(5 / (64 * math.pi)) * (1 - y) ** 2
    assert np.allclose(harmonic, explicit, atol=2e-15)
    amplitude = 1.7 - 0.4j
    recovered = project_mode(amplitude * explicit, projection_row(2, 2, y, weights))
    assert recovered == pytest.approx(amplitude, abs=3e-14)
    leakage = project_mode(explicit, projection_row(3, 2, y, weights))
    assert leakage == pytest.approx(0.0, abs=3e-14)
    harmonic44 = spin_weighted_spherical_harmonic(4, 4, y)
    explicit44 = 3 * math.sqrt(7 / math.pi) * (1 - y) ** 3 * (1 + y) / 16
    assert np.allclose(harmonic44, explicit44, atol=3e-15)
    field = np.vstack((explicit, 2 * explicit))
    assert np.array_equal(scri_slice(field, np.array([0.0, 0.2])), explicit)
    with pytest.raises(OpenBoundaryError):
        scri_slice(field, np.array([0.1, 0.2]))
