
```mermaid
flowchart TB;
    
 A([Sample pool n= 157])--> B[/Rarefy seqs
     @750x50/]   
    subgraph "`**Model A**`"
        B --> C([99 samples])
        C --Subset of 95#--> D[/Randomized samples/]
    end

 A -- Subset #136 --> G[/Randomize samples/]
    subgraph "`**Model B**`"
        G --> H[/Rarefy seqs
     @750x50/]
        H --> I([Unknown #])
    end
    I & D --> T{Decision} --> F[\Downstream analysis/] 

    style D fill:#ffcf11,stroke:#333,stroke-width:4px
    style G fill:#ffcf11,stroke:#333,stroke-width:4px
    style B fill:#ee2d6a,stroke:#333,stroke-width:4px
    style H fill:#ee2d6a,stroke:#333,stroke-width:4px
    style T fill:#aed136,stroke:#333,stroke-width:4px
    style A fill:#31c3e2,stroke:#333,stroke-width:4px
    style C fill:#31c3e2,stroke:#333,stroke-width:4px
    style I fill:#31c3e2,stroke:#333,stroke-width:4px
    linkStyle default stroke: blue
   
