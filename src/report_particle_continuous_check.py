"""Check particle shell sums with a second integration of identical source nodes."""
import argparse
import hashlib
import json
from pathlib import Path
from environment_response import SampledResponse


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    args=parser.parse_args();folder=args.input.parent
    report=json.loads(args.input.read_text());checks=[]
    for shell in report['multipoles']:
        if not shell['complete']:continue
        coherent=0j
        for component in shell['components']:
            path=folder/component['file']
            if hashlib.sha256(path.read_bytes()).hexdigest()!=report['inputs_sha256'][path.name]:
                raise ValueError('Source changed since original particle sum')
            data=json.loads(path.read_text());r0=data['parameters']['metric']['orbital_radius']
            response=SampledResponse.from_report(data)
            radial=response.evaluate([r0])[0,0]
            coherent+=radial*component['equatorial_angular_value']
        old=complex(*shell['particle_field_sum'])
        row=dict(ell=shell['ell'],gauss_field=[old.real,old.imag],
                 continuous_field=[coherent.real,coherent.imag],
                 relative_difference=abs(coherent-old)/abs(old) if old else None)
        checks.append(row);print(row,flush=True)
        result=dict(status='same_source_nodes_not_source_grid_convergence',checks=checks,
            input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
            complete=len(checks)==len(report['complete_shells']))
        output=folder/'particle_field_continuous_check.json'
        temporary=output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(output)


if __name__=='__main__':main()
