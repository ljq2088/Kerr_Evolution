# r0=10M dense-angular source and flux comparison

Running batch: scalar_batch_mg1_08c580cc3167.json. Baseline: scalar_batch_mg1_604b085df0b7.json.
Batch parameters were compared exactly after removing only dense_angular.
The six channels are scalar (ell,m)=(2,2),(4,2),(6,2),(0,0),(2,0),(4,0).
The metric has L18, mg1; angular order18, ordinary radial order8, horizon
logarithmic order64, inner source offset0.0005, source outer320, Green outer1000.
Source grid has 128 nodes. Three workers reconstruct the independent dense
metric cache. Existing default cache and source outputs are retained.

The kappa backend remains the existing finite difference so this run isolates
the angular backend. The C++ radial solver and scalar angular projector are
unchanged. Actual source nodes and weights, source metadata, response fluxes
and complex amplitudes must be checked when the batch completes.

This specifically tests the dominant r10 scalar00 horizon contribution.
A previous r20/mg5 dense comparison cannot establish the effect on this mode.
No flux result or agreement with the paper is claimed yet. Log:
outputs/rp10_mg1_h64_denseangular.log.


The comparison entry point src/report_backend_batch_comparison.py accepts
matching completed default/dense batches, including bound modes with zero
infinity flux. It checks all response parameters after removing only the
angular-backend identifier, identical radial nodes/weights, finite sources,
batch/response flux consistency, and amplitude-squared flux changes.
Its propagating path was regressed against the four completed r20/mg5
channels in dense_mg5_full_flux_comparison.json; relative flux changes agree
within 1e-14. The new incomplete r10 batch is correctly rejected. The bound
path still awaits the completed r10 data; no bound-mode comparison is claimed.


The new report_backend_panel_comparison.py first applies the strict batch
comparison, then differences saved cumulative Green integrals across each
source panel. It checks complete matching panel edges, finite amplitudes,
endpoint agreement with the saved response, and recovery of the total
backend-induced amplitude change. It reports the largest contributing radial
interval and |sum(delta Z_panel)|/sum(|delta Z_panel|). This is attribution
of a discrete difference, not a continuum error or a physical explanation.

The four previously completed r20/mg5 direct channels pass this check in
dense_mg5_radial_panel_comparison.json. The baseline batch was explicitly
restricted to those direct channels (including its conjugate_ells parameter);
the original batch hash and selection are retained. The unmodified full
baseline was correctly rejected because it includes extra conjugate channels.
No r10 panel attribution is claimed until its full batch completes.
