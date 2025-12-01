import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from PROBESt.AI import LogisticRegressionModel, DeepNeuralNetworkModel, TorchClassifier
from PROBESt.filtration import (
    train_filtration_AI, validate_filtration_AI, apply_filtration_AI,
    plot_combined_roc_curves, plot_combined_metrics, plot_learning_curves,
    plot_gail_architecture_comparison
)
from PROBESt.models_registry import (
    ShallowNet, WideNet, ResidualNet, GAILDiscriminator, TabTransformer,
    GAILDeep, GAILWide, GAILNarrow, GAILWithDropout,
    GAILWideDeep, GAILWideDropout, GAILWideBatchNorm, GAILWideExtra, GAILWideBalanced
)

MODELS = {
    "ShallowNet": lambda n: TorchClassifier(ShallowNet(n), weight_pos=5),
    "WideNet": lambda n: TorchClassifier(WideNet(n), weight_pos=5),
    "ResidualNet": lambda n: TorchClassifier(ResidualNet(n), weight_pos=5),
    "GAIL": lambda n: TorchClassifier(GAILDiscriminator(n), weight_pos=5),
    "TabTransformer": lambda n: TorchClassifier(TabTransformer(n), weight_pos=5),
}

def main():
    # Load data
    data_path = 'data/databases/open/test_ML_database.csv'
    data = pd.read_csv(data_path)
    
    # Convert boolean 'type' column to numeric
    data['type'] = data['type'].astype(int)
    
    # Drop character columns but keep the target column
    numeric_columns = data.select_dtypes(include=[np.number]).columns
    data = data[numeric_columns]
    
    # Drop rows with NaN values
    initial_size = len(data)
    data = data.dropna()
    print(f"Dropped {initial_size - len(data)} rows with NaN values")
    
    # Split data into train, validation, and test sets
    train_data, temp_data = train_test_split(data, test_size=0.2, random_state=42)
    val_data, test_data = train_test_split(temp_data, test_size=0.5, random_state=42)
    
    print(f"Training set size: {len(train_data)}")
    print(f"Validation set size: {len(val_data)}")
    print(f"Test set size: {len(test_data)}")
    
    # Get input size
    input_size = train_data.shape[1] - 1
    
    # Train and validate all models from MODELS dictionary
    results = {}
    trained_models = {}
    models_for_plots = {}  # For combined plots: (model, val_data)
    
    for name, constructor in MODELS.items():
        print(f"\n===== Training {name} =====")
        model = constructor(input_size)
        
        # For GAIL, track learning curves during training
        if name == "GAIL":
            X_train = train_data.drop(columns=['type'])
            y_train = train_data['type']
            model.train(X_train, y_train, epochs=100, batch_size=32, 
                       val_data=val_data, track_curves=True)
            trained = model
            # Get training metrics from validation (we'll use val_metrics for results)
            metrics = {}
        else:
            trained, metrics = train_filtration_AI(model, train_data)
        
        val_metrics = validate_filtration_AI(trained, val_data, output_name=f"{name}.png")
        results[name] = val_metrics
        trained_models[name] = trained
        models_for_plots[name] = (trained, val_data)
        
        print(f"\n{name} validation metrics:")
        for metric, value in val_metrics.items():
            print(f"  {metric}: {value:.4f}")
    
    # Train and validate additional baseline models
    print("\n===== Training Deep Neural Network =====")
    dnn_model = DeepNeuralNetworkModel(input_size=input_size)
    trained_dnn_model, dnn_metrics = train_filtration_AI(dnn_model, train_data)
    dnn_val_metrics = validate_filtration_AI(trained_dnn_model, val_data, output_name='DNN.png')
    results["DeepNeuralNetwork"] = dnn_val_metrics
    trained_models["DeepNeuralNetwork"] = trained_dnn_model
    models_for_plots["DeepNeuralNetwork"] = (trained_dnn_model, val_data)
    
    print("\nDeep Neural Network validation metrics:")
    for metric, value in dnn_val_metrics.items():
        print(f"  {metric}: {value:.4f}")
    
    print("\n===== Training Logistic Regression =====")
    lr_model = LogisticRegressionModel()
    trained_lr_model, lr_metrics = train_filtration_AI(lr_model, train_data)
    lr_val_metrics = validate_filtration_AI(trained_lr_model, val_data, output_name='LR.png')
    results["LogisticRegression"] = lr_val_metrics
    trained_models["LogisticRegression"] = trained_lr_model
    models_for_plots["LogisticRegression"] = (trained_lr_model, val_data)
    
    print("\nLogistic Regression validation metrics:")
    for metric, value in lr_val_metrics.items():
        print(f"  {metric}: {value:.4f}")
    
    # Generate combined plots
    print("\n===== Generating combined plots =====")
    output_dir = 'tests_outs'
    plot_combined_roc_curves(models_for_plots, output_dir=output_dir)
    print(f"Combined ROC curves saved to {os.path.join(output_dir, 'combined_roc_curves.png')}")
    
    plot_combined_metrics(results, output_dir=output_dir)
    print(f"Combined metrics plot saved to {os.path.join(output_dir, 'combined_metrics.png')}")
    
    # Select best model based on F1 score
    best_model_name = max(results.keys(), key=lambda m: results[m]["f1"])
    best_model = trained_models[best_model_name]
    
    print(f"\n{'='*60}")
    print(f"Best model: {best_model_name} (F1: {results[best_model_name]['f1']:.4f})")
    print(f"{'='*60}")
    
    # Store GAIL search results for potential use later
    gail_search_results = {}
    gail_trained_models = {}
    
    # Architecture search for GAIL if it's the best model
    if best_model_name == "GAIL":
        print("\n" + "="*60)
        print("GAIL is the best model. Performing architecture search...")
        print("="*60)
        
        # Define GAIL architecture variations
        gail_variations = {
            "GAIL_Deep": lambda n: TorchClassifier(GAILDeep(n), weight_pos=5),
            "GAIL_Wide": lambda n: TorchClassifier(GAILWide(n), weight_pos=5),
            "GAIL_Narrow": lambda n: TorchClassifier(GAILNarrow(n), weight_pos=5),
            "GAIL_Dropout": lambda n: TorchClassifier(GAILWithDropout(n), weight_pos=5),
            "GAIL_Custom1": lambda n: TorchClassifier(GAILDiscriminator(n, hidden1=384, hidden2=192), weight_pos=5),
            "GAIL_Custom2": lambda n: TorchClassifier(GAILDiscriminator(n, hidden1=192, hidden2=96), weight_pos=5),
        }
        
        for variant_name, constructor in gail_variations.items():
            print(f"\n===== Training {variant_name} =====")
            variant_model = constructor(input_size)
            
            # Train with learning curve tracking
            X_train = train_data.drop(columns=['type'])
            y_train = train_data['type']
            variant_model.train(X_train, y_train, epochs=100, batch_size=32, 
                              val_data=val_data, track_curves=True)
            
            # Validate
            variant_val_metrics = validate_filtration_AI(
                variant_model, val_data, output_name=f"{variant_name}.png"
            )
            gail_search_results[variant_name] = variant_val_metrics
            gail_trained_models[variant_name] = variant_model
            
            print(f"\n{variant_name} validation metrics:")
            for metric, value in variant_val_metrics.items():
                print(f"  {metric}: {value:.4f}")
        
        # Compare with original GAIL
        gail_search_results["GAIL_Original"] = results["GAIL"]
        gail_trained_models["GAIL_Original"] = trained_models["GAIL"]
        
        # Find best GAIL variant
        best_gail_variant = max(gail_search_results.keys(), 
                               key=lambda m: gail_search_results[m]["f1"])
        best_gail_model = gail_trained_models[best_gail_variant]
        
        print(f"\n{'='*60}")
        print(f"Best GAIL variant: {best_gail_variant} (F1: {gail_search_results[best_gail_variant]['f1']:.4f})")
        print(f"{'='*60}")
        
        # Plot learning curves for best GAIL variant
        if hasattr(best_gail_model, 'train_losses') and len(best_gail_model.train_losses) > 0:
            print("\nGenerating learning curves for best GAIL variant...")
            plot_learning_curves(best_gail_model, output_dir=output_dir, 
                               output_name=f'learning_curves_{best_gail_variant}.png')
            print(f"Learning curves saved to {os.path.join(output_dir, f'learning_curves_{best_gail_variant}.png')}")
        
        # Update best model if variant is better
        if gail_search_results[best_gail_variant]["f1"] > results["GAIL"]["f1"]:
            print(f"\nBest GAIL variant ({best_gail_variant}) outperforms original GAIL!")
            best_model = best_gail_model
            best_model_name = best_gail_variant
        else:
            print(f"\nOriginal GAIL remains the best.")
        
        # Plot GAIL architecture comparison
        print("\nGenerating GAIL architecture comparison plots...")
        plot_gail_architecture_comparison(gail_search_results, output_dir=output_dir)
        print(f"GAIL comparison plots saved to {os.path.join(output_dir, 'gail_architecture_comparison.png')}")
    
    # Check if best model is any GAIL variant
    is_gail_best = (best_model_name == "GAIL" or 
                    best_model_name.startswith("GAIL") or
                    "GAIL" in best_model_name)
    
    # Further architecture search for any GAIL variant if it's the best model
    if is_gail_best:
        print("\n" + "="*60)
        print(f"{best_model_name} is the best model. Performing further architecture search...")
        print("="*60)
        
        # Determine base architecture type for further search
        if "Wide" in best_model_name:
            # GAIL_Wide variations
            further_variations = {
                "GAIL_Wide_Deep": lambda n: TorchClassifier(GAILWideDeep(n), weight_pos=5),
                "GAIL_Wide_Dropout": lambda n: TorchClassifier(GAILWideDropout(n), weight_pos=5),
                "GAIL_Wide_BatchNorm": lambda n: TorchClassifier(GAILWideBatchNorm(n), weight_pos=5),
                "GAIL_Wide_Extra": lambda n: TorchClassifier(GAILWideExtra(n), weight_pos=5),
                "GAIL_Wide_Balanced": lambda n: TorchClassifier(GAILWideBalanced(n), weight_pos=5),
                "GAIL_Wide_Custom1": lambda n: TorchClassifier(GAILWide(n, hidden1=640, hidden2=320), weight_pos=5),
                "GAIL_Wide_Custom2": lambda n: TorchClassifier(GAILWide(n, hidden1=384, hidden2=192), weight_pos=5),
                "GAIL_Wide_Custom3": lambda n: TorchClassifier(GAILWideDeep(n, hidden1=512, hidden2=256, hidden3=128), weight_pos=5),
            }
        elif "Deep" in best_model_name:
            # GAIL_Deep variations - create deeper and wider versions
            further_variations = {
                "GAIL_Deep_Wide": lambda n: TorchClassifier(GAILDeep(n, hidden1=512, hidden2=256, hidden3=128), weight_pos=5),
                "GAIL_Deep_Dropout": lambda n: TorchClassifier(GAILDeep(n, hidden1=256, hidden2=128, hidden3=64), weight_pos=5),
                "GAIL_Deep_Extra": lambda n: TorchClassifier(GAILDeep(n, hidden1=384, hidden2=192, hidden3=96), weight_pos=5),
                "GAIL_Deep_Custom1": lambda n: TorchClassifier(GAILDeep(n, hidden1=320, hidden2=160, hidden3=80), weight_pos=5),
            }
        elif "Narrow" in best_model_name:
            # GAIL_Narrow variations
            further_variations = {
                "GAIL_Narrow_Deep": lambda n: TorchClassifier(GAILNarrow(n, hidden1=128, hidden2=64), weight_pos=5),
                "GAIL_Narrow_Wide": lambda n: TorchClassifier(GAILNarrow(n, hidden1=256, hidden2=128), weight_pos=5),
                "GAIL_Narrow_Dropout": lambda n: TorchClassifier(GAILNarrow(n, hidden1=128, hidden2=64), weight_pos=5),
            }
        else:
            # Generic GAIL variations for other types
            further_variations = {
                "GAIL_Further_Deep": lambda n: TorchClassifier(GAILDeep(n, hidden1=384, hidden2=192, hidden3=96), weight_pos=5),
                "GAIL_Further_Wide": lambda n: TorchClassifier(GAILWide(n, hidden1=640, hidden2=320), weight_pos=5),
                "GAIL_Further_Dropout": lambda n: TorchClassifier(GAILWithDropout(n, hidden1=384, hidden2=192), weight_pos=5),
                "GAIL_Further_Balanced": lambda n: TorchClassifier(GAILDiscriminator(n, hidden1=256, hidden2=128), weight_pos=5),
            }
        
        further_search_results = {}
        further_trained_models = {}
        
        for variant_name, constructor in further_variations.items():
            print(f"\n===== Training {variant_name} =====")
            variant_model = constructor(input_size)
            
            # Train with learning curve tracking (more epochs for better training)
            X_train = train_data.drop(columns=['type'])
            y_train = train_data['type']
            variant_model.train(X_train, y_train, epochs=150, batch_size=32, 
                              val_data=val_data, track_curves=True)
            
            # Validate
            variant_val_metrics = validate_filtration_AI(
                variant_model, val_data, output_name=f"{variant_name}.png"
            )
            further_search_results[variant_name] = variant_val_metrics
            further_trained_models[variant_name] = variant_model
            
            print(f"\n{variant_name} validation metrics:")
            for metric, value in variant_val_metrics.items():
                print(f"  {metric}: {value:.4f}")
        
        # Get the current best model's metrics
        if best_model_name in gail_search_results:
            original_metrics = gail_search_results[best_model_name]
        else:
            original_metrics = results.get(best_model_name, {})
        
        further_search_results[f"{best_model_name}_Original"] = original_metrics
        further_trained_models[f"{best_model_name}_Original"] = best_model
        
        # Find best variant from further search
        best_further_variant = max(further_search_results.keys(), 
                                  key=lambda m: further_search_results[m].get("f1", 0))
        best_further_model = further_trained_models[best_further_variant]
        
        print(f"\n{'='*60}")
        print(f"Best variant from further search: {best_further_variant} (F1: {further_search_results[best_further_variant].get('f1', 0):.4f})")
        print(f"{'='*60}")
        
        # Combine all GAIL results for comparison plot
        all_gail_results = {}
        if len(gail_search_results) > 0:
            all_gail_results.update(gail_search_results)
        all_gail_results.update(further_search_results)
        
        # Plot GAIL architecture comparison
        print("\nGenerating GAIL architecture comparison plots...")
        plot_gail_architecture_comparison(all_gail_results, output_dir=output_dir)
        print(f"GAIL comparison plots saved to {os.path.join(output_dir, 'gail_architecture_comparison.png')}")
        
        # Plot learning curves for best variant
        if hasattr(best_further_model, 'train_losses') and len(best_further_model.train_losses) > 0:
            print("\nGenerating learning curves for best variant...")
            plot_learning_curves(best_further_model, output_dir=output_dir, 
                               output_name=f'learning_curves_{best_further_variant}.png')
            print(f"Learning curves saved to {os.path.join(output_dir, f'learning_curves_{best_further_variant}.png')}")
        
        # Update best model if variant is better
        current_f1 = further_search_results.get(f"{best_model_name}_Original", {}).get("f1", 0)
        if further_search_results[best_further_variant].get("f1", 0) > current_f1:
            print(f"\nBest variant ({best_further_variant}) outperforms original!")
            best_model = best_further_model
            best_model_name = best_further_variant
        else:
            print(f"\nOriginal {best_model_name} remains the best.")
        
        # Train the best GAIL model for extended epochs
        print(f"\n{'='*60}")
        print(f"Training {best_model_name} for extended epochs (300 epochs)...")
        print(f"{'='*60}")
        X_train = train_data.drop(columns=['type'])
        y_train = train_data['type']
        best_model.train(X_train, y_train, epochs=300, batch_size=32, 
                        val_data=val_data, track_curves=True)
        
        # Re-validate after extended training
        final_val_metrics = validate_filtration_AI(
            best_model, val_data, output_name=f"{best_model_name}_final.png"
        )
        print(f"\nFinal {best_model_name} validation metrics after extended training:")
        for metric, value in final_val_metrics.items():
            print(f"  {metric}: {value:.4f}")
        
        # Plot final learning curves
        if hasattr(best_model, 'train_losses') and len(best_model.train_losses) > 0:
            print("\nGenerating final learning curves...")
            plot_learning_curves(best_model, output_dir=output_dir, 
                               output_name=f'learning_curves_{best_model_name}_final.png')
            print(f"Final learning curves saved to {os.path.join(output_dir, f'learning_curves_{best_model_name}_final.png')}")
    
    # Plot learning curves for the best model (if it has learning curve data)
    if hasattr(best_model, 'train_losses') and len(best_model.train_losses) > 0:
        print("\nGenerating learning curves for best model...")
        plot_learning_curves(best_model, output_dir=output_dir, 
                           output_name='learning_curves_best_model.png')
        print(f"Learning curves saved to {os.path.join(output_dir, 'learning_curves_best_model.png')}")
    
    # Apply best model to test set
    print("\nApplying best model to test set...")
    test_predictions = apply_filtration_AI(best_model, test_data)
    
    # Save predictions
    os.makedirs(output_dir, exist_ok=True)
    test_predictions.to_csv(os.path.join(output_dir, 'test_predictions.csv'), index=False)
    print(f"\nPredictions saved to {os.path.join(output_dir, 'test_predictions.csv')}")

if __name__ == '__main__':
    main() 