declare namespace $3Dmol {
  function createViewer(element: HTMLElement, config?: any): Viewer;
  
  interface Viewer {
    addModel(data: string, format: string): Model;
    setStyle(sel: any, style: any): void;
    zoomTo(): void;
    zoom(factor: number): void;
    spin(enable: boolean | number, axis?: any): void;
    render(): void;
    addLine(options: any): void;
    addArrow(options: any): void;
    addLabel(text: string, options: any): void;
    clear(): void;
    resize(): void;
    getModel(): Model;
  }
  
  interface Model {
    setStyle(sel: any, style: any): void;
    setColorByProperty(sel: any, prop: string, scheme: any): void;
    addPolyhedra(options: any): void;
  }
  
  namespace Gradient {
    class RWB {
      constructor(min: number, max: number);
    }
    class ROYGB {
      constructor(min: number, max: number);
    }
  }
}

declare interface Window {
  $3Dmol: typeof $3Dmol;
} 