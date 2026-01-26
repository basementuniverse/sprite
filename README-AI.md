# Sprite Library API Reference (AI Documentation)

## Package
`@basementuniverse/sprite` - 2D sprite component with animations, directions, and attachment points

## Core Types

### vec2
```typescript
{ x: number; y: number }
```

### SpriteOptions
```typescript
{
  position?: vec2;                    // Default: {x:0, y:0}
  size?: vec2;                        // Default: base image size or first animation frame size or {x:0, y:0}
  origin?: vec2;                      // Default: center of sprite (size * 0.5)
  scale?: number;                     // Default: 1
  rotation?: number;                  // Radians, default: 0
  directions: string[];               // Array of valid direction names, default: ['default']
  defaultDirection: string;           // Default: 'default'
  image?: HTMLImageElement | HTMLCanvasElement; // Fallback base image
  animations: {
    [name: string]: {
      [direction: string]: SpriteAnimationOptions;
    };
  };
  defaultAnimation: string;
  attachmentPoints?: SpriteAttachmentPointOptions[];
  preRender?: (context: CanvasRenderingContext2D, sprite: Sprite) => void;
  postRender?: (context: CanvasRenderingContext2D, sprite: Sprite) => void;
  debug?: Partial<SpriteDebugOptions> | boolean;
}
```

### SpriteAnimationOptions
```typescript
{
  name: string;
  frameCount?: number;                // Default: 1
  frameRate?: number;                 // FPS, default: 1
  mode?: SpriteAnimationRepeatMode;   // Default: Repeat
  images?: (HTMLImageElement | HTMLCanvasElement)[];
  attachmentPointKeyframes?: {
    [attachmentPointName: string]: SpriteAttachmentPointKeyframe[];
  };
}
```

### SpriteAnimationRepeatMode (enum)
```typescript
'repeat'              // Loop indefinitely
'play-once-and-stop'  // Play once, stop on last frame
'play-once-and-reset' // Play once, reset to first frame
```

### SpriteAttachmentPointOptions
```typescript
{
  name: string;
  offset: vec2;  // Position offset from sprite origin
}
```

### SpriteAttachmentPointKeyframe
```typescript
{
  frame: number;  // 0-based frame index
  offset: vec2;   // Attachment point offset for this keyframe
}
```

### SpriteDebugOptions
```typescript
{
  showSpriteTransforms: boolean;
  showSpriteBoundingBox: boolean;
  showAttachmentPoints: boolean;
}
```

## Sprite Class

### Constructor
```typescript
new Sprite(options?: Partial<SpriteOptions>)
```

### Public Properties
```typescript
sprite.position: vec2      // Read/write
sprite.size: vec2          // Read/write
sprite.origin: vec2        // Read/write
sprite.scale: number       // Read/write
sprite.rotation: number    // Read/write (radians)
sprite.direction: string   // Read/write getter/setter
sprite.animation: string   // Read/write getter/setter
```

### Public Methods
```typescript
sprite.update(dt: number): void
// Updates animation state, current frame, current image, attachment points

sprite.draw(context: CanvasRenderingContext2D): void
// Renders sprite with transforms, calls preRender/postRender hooks, draws debug info

sprite.playAnimation(): void
// Resumes animation playback

sprite.pauseAnimation(): void
// Pauses animation playback

sprite.resetAnimation(): void
// Resets animation to frame 0

sprite.getAttachmentPoint(name: string): vec2 | null
// Returns current world-space position of named attachment point
```

## Content Manager Integration

### SpriteOptionsData
```typescript
Partial<Omit<SpriteOptions, 'image' | 'preRender' | 'postRender' | 'debug'>> & {
  imageName?: string;  // Instead of 'image', fetched from ContentManager
  animations?: {
    [name: string]: {
      [direction: string]: SpriteAnimationOptionsData;
    };
  };
}
```

### SpriteAnimationOptionsData
```typescript
Omit<SpriteAnimationOptions, 'images'> & {
  imageNames?: string[];  // Instead of 'images', fetched from ContentManager
}
```

### Content Processor
```typescript
async function spriteOptionsContentProcessor(
  content: Record<string, {name: string; type: string; content: any; status: string}>,
  data: {name: string; type: string; content: any; status: string}
): Promise<void>
```
Converts SpriteOptionsData to SpriteOptions by loading images from ContentManager

## Type Guards
```typescript
function isSpriteOptionsData(value: unknown): value is SpriteOptionsData
```

## Animation Direction Fallback
- Direction '*' in animations means "available for all directions"
- Lookup order: specific direction → '*' → first available direction

## Attachment Points
- Defined with default offset from sprite origin
- Can be animated per-frame via keyframes in each animation
- Values interpolated between keyframes (ascending frame order)
- Retrieved via `getAttachmentPoint(name)` returns current world-space position

## Rendering Behavior
- Sprite transforms: translate(position) → scale(scale) → rotate(rotation)
- Image drawn at (-origin.x, -origin.y) relative to transformed origin
- If animation frame image missing: falls back to base image → null (no render)
- Attachment points and debug markers drawn in transformed space

## Animation State
- `playing` boolean controls whether animation advances
- `currentFrame` (0-based index) and `currentFrameTime` track progress
- Frame advances when currentFrameTime >= (1 / frameRate)
- Repeat modes control behavior at end: loop, stop on last, or reset to first
- Switching animations clamps currentFrame if new animation has fewer frames

## Default Values
```typescript
// Default sprite has:
directions: ['default']
defaultDirection: 'default'
defaultAnimation: 'default'
animations: {
  default: {
    '*': {
      name: 'default',
      frameCount: 1,
      frameRate: 1,
      mode: 'play-once-and-stop'
    }
  }
}
```
