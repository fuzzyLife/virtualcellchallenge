"""Apply a saved neural model using exactly the named 300 measured genes."""
import argparse
from pathlib import Path
from neural_runtime import activate_local_runtime
activate_local_runtime()
import torch
from threadpoolctl import threadpool_limits
from predict_gene300_models import read_inputs
from train_gene300_neural import OUT, load_network, write_output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--context', required=True, choices=['A', 'B', 'C'])
    parser.add_argument('--objective', required=True, choices=['normalized', 'raw'])
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--model-dir', type=Path, default=OUT)
    parser.add_argument('--device', choices=['cuda', 'cpu'], default='cpu')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f'Refusing to overwrite {args.output}')
    model, network = load_network(args.model_dir / f'model_{args.objective}_{args.context}.npz', args.device)
    if str(model['context']) != args.context or str(model['objective']) != args.objective:
        raise ValueError('Model context/objective mismatch')
    inputs, obs = read_inputs(args.input, model, args.context)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(8)
    with threadpool_limits(limits=8):
        write_output(args.output, model, network, inputs, obs, 'new_input_all_controls_model')
    print(f'Saved {len(inputs):,} cells x {len(model["output_genes"]):,} genes; objective={args.objective}; output={args.output}')
